% heatmaps to show tuning to 'micro' structure of task: steps, licks, wisk contacts

unitInfo = getUnitInfo('nucleiOnly', true);
unitInfo = unitInfo(1:20,:);  % !!! temp
sessions = unique(unitInfo.session);
paper2_config;


%% compute neural and behavioral responses for steps, licks, and whisker contactsd
% this is slow because requires interpolating over variable duraction steps
% can load E:\lab_files\paper2\intermediate_analysis\microheatmaps.mat if
% already computed...

% settings
x.step = linspace(0, 1, 100);
x.lick = linspace(-.15, .15, 100);
x.wisk = linspace(-.1, .4, 100);
maxStepDuration = .5;


% create unitInfo and sesInfo tables with initialized variables
sesInfo = table('RowNames', sessions);
for dv = {'step', 'lick', 'wisk'}
    if strcmp(dv, 'step')
        sesInfo.([dv{1} '_mean']) = nan(length(sessions), 3, length(x.(dv{1})));  % x y z trajectories for each step
    else
        sesInfo.([dv{1} '_mean']) = nan(length(sessions), length(x.(dv{1})));  % x y z trajectories for each step
    end
    unitInfo.([dv{1} '_resp']) = nan(height(unitInfo), length(x.(dv{1})));  % x y z trajectories for each step
    unitInfo.([dv{1} '_respFull']) = cell(height(unitInfo), 1);  % x y z trajectories for each step
    % todo: save full response as mean, not just mean per unit?
end
unitInfo.mean = nan(height(unitInfo), 1);
unitInfo.std = nan(height(unitInfo), 1);



for i = 1:height(sesInfo)
    
    disp(i/length(sessions))
    
    % load session Data
    p = load(['E:\lab_files\paper2\modelling\predictors\' sessions{i} '_predictors.mat']);
    e = load(['E:\lab_files\paper2\modelling\neuralData\' sessions{i} '_neuralData.mat']);
    
    xyz = p.predictors{{'paw4RH_x', 'paw4RH_y', 'paw4RH_z'}, 'data'}; xyz = cat(1, xyz{:});
    wisk = p.predictors{'whiskerAngle', 'data'}{1};
    jaw = p.predictors{'jaw', 'data'}{1};
    t = p.predictors{'jaw', 't'}{1};  % assumes all predictors share time axis (should be true unless something upstream changes)!
    spkRates = interp1(e.timeStamps, e.spkRates', t)'; % put fr on same time axis as lick rate and vel
    
    
    events.wisk = p.predictors{'whiskerContact', 'data'}{1};
    events.lick = p.predictors{'lick', 'data'}{1};
    events.step = p.predictors{'paw4RH_stride', 'data'}{1};
    validSteps = diff(events.step, [], 2) < maxStepDuration;  % remove steps that are too long
    events.step = events.step(validSteps, :);
    
    % compute neural and behavioral responses for...
    
    % wisk
    responses.wisk = permute(interp1(t, [spkRates; wisk]', events.wisk + x.wisk), [3,2,1]);
    responses.lick = permute(interp1(t, [spkRates; jaw]', events.lick + x.lick), [3,2,1]);
    
    stacked = [spkRates; xyz];
    responses.step = nan(size(stacked,1), length(x.step), size(events.step,1));
    for j = 1:size(events.step,1)  % for each step
        bins = (t>=events.step(j,1)-.5) & (t<=events.step(j,2)+.5);  % +-.5 buffer to include edge samples in interpolation :)
        stepx = linspace(events.step(j,1), events.step(j,2), length(x.step));
        responses.step(:,:,j) = interp1(t(bins), stacked(:,bins)', stepx)';
    end
    
    % save behavioral responses in sesInfo
    for dv = {'step', 'lick', 'wisk'}
        if strcmp(dv{1}, 'step')
            sesInfo{i, [dv{1} '_mean']} = permute(nanmean(responses.(dv{1})(end-2:end, :, :), 3), [3 1 2]);  % the permute is a trick to add a leading singleton dimension
        else
            sesInfo{i, [dv{1} '_mean']} = nanmean(responses.(dv{1})(end, :, :), 3);
        end
    end
    
    % save neural responses in unitInfo
    unitInds = find(strcmp(unitInfo.session, sessions{i}) & ismember(unitInfo.unit, e.unit_ids));
    for j = 1:length(unitInds)
        sesInd = find(e.unit_ids == unitInfo.unit(unitInds(j)));  % index for unit within session
        for dv = {'step', 'lick', 'wisk'}
            unitInfo{unitInds(j), [dv{1} '_respFull']} = {squeeze(responses.(dv{1})(sesInd,:,:))'};  % matrix of response across trials
            unitInfo{unitInds(j), [dv{1} '_resp']} = nanmean(responses.(dv{1})(sesInd,:,:), 3);      % mean response
            
            % unit mean and std (for z scoring later)
            unitInfo.mean(unitInds(j)) = nanmean(e.spkRates(sesInd,:));
            unitInfo.std(unitInds(j)) = nanstd(e.spkRates(sesInd,:));
        end
    end
end

save('E:\lab_files\paper2\intermediate_analysis\microheatmaps.mat', 'sesInfo', 'unitInfo')

%% get response significance

% settings
wiskPre  = [-.15 .05];  % window for pre and post periods for wisk stats
wiskPost = [0 .1];
p_thresh = .01 / height(unitInfo);  % bonferroni corrected p value

% inits
wiskPreBins  = x.wisk>=wiskPre(1)  & x.wisk<=wiskPre(2);
wiskPostBins = x.wisk>=wiskPost(1) & x.wisk<=wiskPost(2);
unitInfo.wisk_sig = nan(height(unitInfo), 1);

for i = 1:height(unitInfo)
    
    % wisk
    resp = unitInfo.wisk_respFull{i};
    pre  = nanmean(resp(:, wiskPreBins),2);
    post = nanmean(resp(:, wiskPostBins),2);
    unitInfo.wisk_sig(i) = ranksum(pre, post);
end


%% plot :)

% todo: limit to significant responses // sort to look dope!

% settings
% lims = [-1 1];  % (z-score)

close all; figure('color', 'white', 'menubar', 'none', 'position', [687.00 2.00 382.00 1408.00])

% steps
subplot(9,1,1); hold on
mask = linspace(.5, 1, 3);
lns = nan(1,3);
for j = 1:3
    resp = squeeze(sesInfo.step_mean(:,j,:));
    c = cfg.velColor * mask(j);
    errorplot(x.step, resp - nanmean(resp,2), 'color', c, 'linewidth', 3)
    lns(j) = plot(nan, nan, 'Color', c, 'LineWidth', 3);
end
legend(lns, {'x', 'y', 'z'}, 'AutoUpdate', 'off');
ylabel('paw position');

addylabel('paw position')
addscalebar(.01, 'y', '10 mm');
set(gca, 'xlim', [x.step(1) x.step(end)], 'xcolor', 'none', 'ycolor', 'none')

subplot(9,1,2:3); hold on
resp = (unitInfo.step_resp - unitInfo.mean) ./ unitInfo.std;  % z score based on overall fr
clims = max(abs(resp(:))) * [-1 1];
[~, maxInds] = max(resp, [], 2);
[~, sortInds] = sort(maxInds, 'descend');
imagesc(x.step, 1:size(resp,1), resp(sortInds,:), clims)
colormap(cfg.heatmapColors)
xlabel('fraction of step');
set(gca, 'XLim', [x.step(1) x.step(end)], 'ylim', [.5 size(resp,1)+.5], 'ycolor', 'none', cfg.axArgs{:}); limitticks
addylabel('firing rate (z-score)')

% licks
subplot(9,1,4); hold on
errorplot(x.lick, sesInfo.lick_mean, 'color', cfg.lickColor, 'linewidth', 3)
set(gca, 'xcolor', 'none', 'ycolor', 'none', 'xlim', [x.lick(1) x.lick(end)])
addylabel('jaw position')
addscalebar(.001, 'y', '1 mm');

subplot(9,1,5:6); hold on
resp = (unitInfo.lick_resp - unitInfo.mean) ./ unitInfo.std;
clims = max(abs(resp(:))) * [-1 1];
[~, maxInds] = max(resp, [], 2);
[~, sortInds] = sort(maxInds, 'descend');
imagesc(x.lick, 1:size(resp,1), resp(sortInds,:), clims)
set(gca, 'XLim', [x.lick(1) x.lick(end)], 'ylim', [.5 size(resp,1)+.5], 'ycolor', 'none', cfg.axArgs{:}); limitticks
xlabel('time from lick (s)');
addylabel('firing rate (z-score)')

% wisk
subplot(9,1,7); hold on
errorplot(x.wisk, sesInfo.wisk_mean, 'color', cfg.wiskColor, 'linewidth', 3)
set(gca, 'xcolor', 'none', 'ycolor', 'none', 'xlim', [x.wisk(1) x.wisk(end)])
addylabel('whisker angle')
addscalebar(10, 'y', '10 \circ');

subplot(9,1,8:9); hold on
resp = (unitInfo.wisk_resp - unitInfo.mean) ./ unitInfo.std;
resp = resp(unitInfo.wisk_sig < p_thresh, :);
clims = max(abs(resp(:))) * [-1 1];

% colbins = x.wisk>-.05 & x.wisk<.05;
colbins = x.wisk>-0 & x.wisk<.1;
[~, sortInds] = sort(nanmean(resp(:, colbins),2));

% groups = clusterResponses(resp, 'plot', false, 'nclusters', 4, 'pcs', 2);
% respMag = resp(:, find(x.wisk>.06,1,'first'));
% [~, sortInds] = sortrows([groups respMag], 'ascend');

imagesc(x.wisk, 1:size(resp,1), resp(sortInds,:), clims)
set(gca, 'XLim', [x.wisk(1) x.wisk(end)], 'ylim', [.5 size(resp,1)+.5], 'ycolor', 'none', cfg.axArgs{:}); limitticks
xlabel('time from whisker contact (s)');
addylabel('firing rate (z-score)')
























