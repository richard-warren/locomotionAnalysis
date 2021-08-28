% heatmaps to show tuning to 'micro' structure of task: steps, licks, wisk contacts

unitInfo = getUnitInfo('nucleiOnly', true, 'frstats', true);
% unitInfo = unitInfo(1:20,:);  % !!! temp
sessions = unique(unitInfo.session);
paper2_config;


%% settings
x.step = linspace(0, 1, 100);
x.lick = linspace(-.15, .15, 100);
x.wisk = linspace(-.1, .4, 100);
maxStepDuration = .5;
maxLickInterval = .5;

wiskPre  = [-.15 .05];  % window for pre and post periods for wisk stats
wiskPost = [0 .1];
p_thresh = .05 / height(unitInfo);  % bonferroni corrected p value


%% compute neural and behavioral responses for steps, licks, and whisker contactsd
% this is slow because requires interpolating over variable duration steps
% can load E:\lab_files\paper2\intermediate_analysis\microheatmaps.mat if
% already computed...




% inits
wiskPreBins  = x.wisk>=wiskPre(1)  & x.wisk<=wiskPre(2);
wiskPostBins = x.wisk>=wiskPost(1) & x.wisk<=wiskPost(2);

% create unitInfo and sesInfo tables with initialized variables
sesInfo = table('RowNames', sessions);
for dv = {'step', 'lick', 'wisk'}
    if strcmp(dv, 'step')
        sesInfo.([dv{1} '_mean']) = nan(length(sessions), 3, length(x.(dv{1})));
    else
        sesInfo.([dv{1} '_mean']) = nan(length(sessions), length(x.(dv{1})));
    end
    unitInfo.([dv{1} '_resp']) = nan(height(unitInfo), length(x.(dv{1})));
    unitInfo.([dv{1} '_respFull']) = cell(height(unitInfo), 1);
    unitInfo.([dv{1} '_sig']) = nan(height(unitInfo), 1);
    unitInfo.([dv{1} '_k']) = nan(height(unitInfo), 1);
    unitInfo.([dv{1} '_phase']) = cell(height(unitInfo), 1);  % phase at which spikes occurred
end



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
    
    % get event times
    events.wisk = p.predictors{'whiskerContact', 'data'}{1};
    
    events.lick = p.predictors{'lick', 'data'}{1};
    events.lickEpoch = [events.lick(1:end-1) events.lick(2:end)];  % turn into epoch
    validLicks = diff(events.lickEpoch, [], 2) < maxLickInterval;  % remove steps that are too long
    events.lickEpoch = events.lickEpoch(validLicks, :);
    
    events.step = p.predictors{'paw4RH_stride', 'data'}{1};
    validSteps = diff(events.step, [], 2) < maxStepDuration;  % remove steps that are too long
    events.step = events.step(validSteps, :);
    
    % compute neural and behavioral responses for...
    
    % wisk
    responses.wisk = permute(interp1(t, [spkRates; wisk]', events.wisk + x.wisk), [3,2,1]);
    
    % lick
    responses.lick = permute(interp1(t, [spkRates; jaw]', events.lick + x.lick), [3,2,1]);
    
    phase.lick = nan(size(t));
    for j = 1:size(events.lickEpoch,1)  % for each lick
        bins = (t>=events.lickEpoch(j,1)) & (t<events.lickEpoch(j,2));
        phase.lick(bins) = interp1(events.lickEpoch(j,:), [0 2*pi], t(bins), 'linear', 'extrap');
    end
    
    % step
    stacked = [spkRates; xyz];
    phase.step = nan(size(t));
    responses.step = nan(size(stacked,1), length(x.step), size(events.step,1));
    for j = 1:size(events.step,1)  % for each step
        bins = (t>=events.step(j,1)-.1) & (t<=events.step(j,2)+.1);  % +-.1 buffer to include edge samples in interpolation :)
        stepx = linspace(events.step(j,1), events.step(j,2), length(x.step));
        responses.step(:,:,j) = interp1(t(bins), stacked(:,bins)', stepx)';
        
        bins = (t>=events.step(j,1)) & (t<events.step(j,2));  % +-.5 buffer to include edge samples in interpolation :)
        phase.step(bins) = interp1(events.step(j,:), [0 2*pi], t(bins), 'linear', 'extrap');
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
            
            % statistical tests
            if strcmp(dv{1}, 'wisk')
                resp = unitInfo.wisk_respFull{unitInds(j)};
                pre  = nanmean(resp(:, wiskPreBins),2);
                post = nanmean(resp(:, wiskPostBins),2);
                unitInfo{unitInds(j), [dv{1} '_sig']} = ranksum(pre, post);
            else
                phaseAtSpks = interp1(t, phase.(dv{1}), e.spkTimes{sesInd}, 'nearest');
                phaseAtSpks = phaseAtSpks(~isnan(phaseAtSpks));
                [unitInfo{unitInds(j), [dv{1} '_k']}, unitInfo{unitInds(j), [dv{1} '_sig']}] = ...
                    ktestnull(phaseAtSpks, 'res', 40, 'bootstraps', 100, 'plot', false);
                unitInfo{unitInds(j), [dv{1} '_phase']} = {phaseAtSpks};  % save phase at which spikes occured
            end
        end
    end
end


% modulation index
for dv = {'step', 'lick'}
    maxs = max(unitInfo.([dv{1} '_resp']), [], 2);
    mins = min(unitInfo.([dv{1} '_resp']), [], 2);
    unitInfo.([dv{1} '_mi']) = (maxs-mins) ./ (maxs + mins);
end

save('E:\lab_files\paper2\intermediate_analysis\microheatmaps.mat', 'sesInfo', 'unitInfo', '-v7.3')





%% plot :)

% settings
% k_thresh = .05;


% close all; figure('color', 'white', 'menubar', 'none', 'position', [846.00 92.00 248.00 1079.00])
close all; figure('color', 'white', 'menubar', 'none', 'position', [243.00 567.00 724.00 608.00])
nrows = nan(1,3);  % number of rows for each plot

% steps
subplot(3,3,1); hold on
% subplot(9,1,1); hold on
mask = linspace(.5, 1, 3);  % creates color gradient from cfg.velColor
lns = nan(1,3);
for j = 1:3
    resp = squeeze(sesInfo.step_mean(:,j,:));
    c = cfg.velColor * mask(j);
    errorplot(x.step, resp - nanmean(resp,2), 'color', c, 'linewidth', 2)
    lns(j) = plot(nan, nan, 'Color', c, 'LineWidth', 2);
end
legend(lns, {'x', 'y', 'z'}, 'AutoUpdate', 'off', 'FontName', cfg.font, 'FontSize', cfg.fontsize);
addylabel('paw position', cfg.font, cfg.fontsize)
addscalebar(.01, 'y', '10 mm', cfg.font, cfg.fontsize);
set(gca, 'xlim', [x.step(1) x.step(end)], 'xcolor', 'none', 'ycolor', 'none', cfg.axArgs{:})

subplot(3,3,[4 7]); hold on
% subplot(9,1,2:3); hold on

mn = cellfun(@(x) nanmean(x(:)), unitInfo.step_respFull);
st = cellfun(@(x) nanstd(x(:)), unitInfo.step_respFull);
resp = (unitInfo.step_resp - mn) ./ st;  % z score based on response period firing rate
% resp = (unitInfo.step_resp - unitInfo.mean) ./ unitInfo.std;  % z score based on overall fr
resp = resp(unitInfo.step_sig<p_thresh, :);  % sig limit
nrows(1) = size(resp,1);
clims = max(abs(resp(:))) * [-1 1];
[~, maxInds] = max(resp, [], 2); [~, sortInds] = sort(maxInds, 'descend');

imagesc(x.step, 1:size(resp,1), resp(sortInds,:), clims)
colormap(cfg.heatmapColors)
xlabel('fraction of step');
set(gca, 'XLim', [x.step(1) x.step(end)], 'ylim', [.5 size(resp,1)+.5], 'ycolor', 'none', cfg.axArgs{:}); limitticks
addylabel('firing rate (z-score)', cfg.font, cfg.fontsize)

% licks
subplot(3,3,2); hold on
% subplot(9,1,4); hold on
errorplot(x.lick, sesInfo.lick_mean*1000, 'color', cfg.lickColor, 'linewidth', 2)  % for some fucked up reason multiplying by 1000 causes the text to be editable when exporting to illustrating... took me 2 hours to figure this our... :(
set(gca, 'xcolor', 'none', 'ycolor', 'none', 'xlim', [x.lick(1) x.lick(end)], cfg.axArgs{:})
addylabel('jaw position', cfg.font, cfg.fontsize)
addscalebar(.001*1000, 'y', '1 mm', cfg.font, cfg.fontsize);

subplot(3,3,[5 8]); hold on
% subplot(9,1,5:6); hold on
mn = cellfun(@(x) nanmean(x(:)), unitInfo.lick_respFull);
st = cellfun(@(x) nanstd(x(:)), unitInfo.lick_respFull);
resp = (unitInfo.lick_resp - mn) ./ st;  % z score based on response period firing rate
% resp = (unitInfo.lick_resp - unitInfo.mean) ./ unitInfo.std;
resp = resp(unitInfo.lick_sig<p_thresh, :);  % sig limit
nrows(2) = size(resp,1);
clims = max(abs(resp(:))) * [-1 1];
[~, sortInds] = sort(resp(:, find(x.lick>0,1,'first')));
imagesc(x.lick, 1:size(resp,1), resp(sortInds,:), clims)
xlabel('time from lick (s)');
set(gca, 'XLim', [x.lick(1) x.lick(end)], 'xtick', [x.lick(1) 0 x.lick(end)], ...
    'ylim', [.5 size(resp,1)+.5], 'ycolor', 'none', cfg.axArgs{:});
% addylabel('firing rate (z-score)', cfg.font, cfg.fontsize)

% wisk
subplot(3,3,3); hold on
% subplot(9,1,7); hold on
errorplot(x.wisk, sesInfo.wisk_mean, 'color', cfg.wiskColor, 'linewidth', 2);
set(gca, 'xcolor', 'none', 'ycolor', 'none', 'xlim', [x.wisk(1) x.wisk(end)])
addylabel('whisker angle', cfg.font, cfg.fontsize);
addscalebar(10, 'y', '10 degrees', cfg.font, cfg.fontsize);

subplot(3,3,[6 9]); hold on
% subplot(9,1,8:9); hold on
mn = cellfun(@(x) nanmean(x(:)), unitInfo.wisk_respFull);
st = cellfun(@(x) nanstd(x(:)), unitInfo.wisk_respFull);
resp = (unitInfo.wisk_resp - mn) ./ st;  % z score based on response period firing rate
% resp = (unitInfo.wisk_resp - unitInfo.mean) ./ unitInfo.std;
resp = resp(unitInfo.wisk_sig < p_thresh, :);
nrows(3) = size(resp,1);
clims = max(abs(resp(:))) * [-1 1];

colbins = x.wisk>-0 & x.wisk<.1;
[~, sortInds] = sort(nanmean(resp(:, colbins),2));

imagesc(x.wisk, 1:size(resp,1), resp(sortInds,:), clims);
set(gca, 'XLim', [x.wisk(1) x.wisk(end)], 'xtick', [0 x.wisk(end)], ...
    'ylim', [.5 size(resp,1)+.5], 'ycolor', 'none', cfg.axArgs{:});
xlabel({'time from', 'whisker contact (s)'});
% addylabel('firing rate (z-score)', cfg.font, cfg.fontsize)



% adjust heatmaps to reflect number of rows
pos = get(gca, 'position');
hgts = (nrows / max(nrows)) * pos(4);  % target heights for each heatmap
for i = 1:3
    subplot(3,3,[3 6]+i);
    pos = get(gca, 'position');
    set(gca, 'Position', [pos(1) pos(2)+(max(hgts-hgts(i))) pos(3) hgts(i)])
end

% adjust heatmaps to reflect number of rows
% pos = get(gca, 'position');
% hgts = (nrows / max(nrows)) * pos(4);  % target heights for each heatmap
% for i = 1:3
%     subplot(9,1,[-1 0]+i*3);
%     pos = get(gca, 'position');
%     set(gca, 'Position', [pos(1) pos(2)+(max(hgts)-hgts(i)) pos(3) hgts(i)])
% end


set(gcf, 'Renderer', 'painters')
saveas(gcf, 'E:\lab_files\paper2\paper_figures\matlab\microheatmaps', 'svg')



% %% check circular stats against PSTHs
% 
% rows = 1:5;
% 
% close all
% figure('color', 'white', 'position', [2.00 2.00 1278.00 1354.00])
% % ks = nan(1, length(rows));
% 
% for i = rows
%     
%     % fr psth
%     subplot(length(rows), 3, (i-1)*3+1)
%     errorplot(x.step, unitInfo{i, 'step_respFull'}{1})
%     
%     % phase histo
%     subplot(length(rows), 3, (i-1)*3+2)
%     p = unitInfo.step_phase{i};
%     histogram(p, 100, 'Normalization', 'probability');
%     set(gca, 'xlim', [0 2*pi])
%     
%     % kuipers
%     subplot(length(rows), 3, (i-1)*3+3)
% %     [sig, k] = circ_kuipertest(linspace(0,2*pi,500), p, 40, true);
%     k = kuipertest2(linspace(0,2*pi,500), p, 40, true);
% end



%% example units

% settings
maxtrials = 200;
scatsz = 2;

units.step.session = {'200818_000', '200201_000', '201208_000'};
units.step.unit    = [45 88 191];
units.lick.session = {'200708_000', '200817_000', '200902_000'};
units.lick.unit    = [11 223 206];

units.wisk.session = {'200116_000', '201102_000', '201208_000'};
units.wisk.unit    = [30 254 199];



% inits
close all
figure('color', 'white', 'position', [523.00 592.00 543.00 466.00], 'menubar', 'none')
nunits = length(units.step.session);  % assumes same numer of units for all dvs!

% step
for i = 1:nunits
    session = units.step.session{i};
    unit = units.step.unit(i);
    
    e = load(['E:\lab_files\paper2\modelling\neuralData\' session '_neuralData.mat']);
    spktimes = e.spkTimes{e.unit_ids==unit};
    spkrate  = e.spkRates(e.unit_ids==unit,:);
    tstart = e.timeStamps(find(~isnan(spkrate), 1, 'first'));
    tend = e.timeStamps(find(~isnan(spkrate), 1, 'last'));
    
    p = load(['E:\lab_files\paper2\modelling\predictors\' session '_predictors.mat']);
    evts = p.predictors{'paw4RH_stride', 'data'}{1};
    validSteps = diff(evts, [], 2) < maxStepDuration;  % remove steps that are too long
    evts = evts(validSteps, :);
    evts = evts(evts(:,1)>tstart & evts(:,2)<tend, :);
    
%     subplots = [6,ncols,i+0*ncols; 6,ncols,i+1*ncols];
    subplots = [nunits*2,3,(i-1)*6+1; nunits*2,3,(i-1)*6+4];
    plotPSTH(spktimes, evts, 'kernel', 'gauss', 'eventLims', [x.lick(1) x.lick(end)], ...
        'color', cfg.velColor, 'subplots', subplots, 'maxEpochs', maxtrials, 'scatSize', scatsz, 'binNum', 100); rasterize
    if i==nunits; xlabel('fraction of step cycle'); else; set(gca, 'XTickLabel', []); end
    set(gca, 'xtick', [x.step(1) x.step(end)])
end

% lick
for i = 1:nunits
    session = units.lick.session{i};
    unit = units.lick.unit(i);
    
    e = load(['E:\lab_files\paper2\modelling\neuralData\' session '_neuralData.mat']);
    spktimes = e.spkTimes{e.unit_ids==unit};
    spkrate  = e.spkRates(e.unit_ids==unit,:);
    tstart = e.timeStamps(find(~isnan(spkrate), 1, 'first'));
    tend = e.timeStamps(find(~isnan(spkrate), 1, 'last'));
    
    p = load(['E:\lab_files\paper2\modelling\predictors\' session '_predictors.mat']);
    evts = p.predictors{'lick', 'data'}{1};
    evts = evts(evts>tstart & evts<tend);
    evts = evts(sort(randperm(length(evts), min(maxtrials, length(evts)))), :);
    
%     subplots = [6,ncols,i+2*ncols; 6,ncols,i+3*ncols];
    subplots = [nunits*2,3,(i-1)*6+2; nunits*2,3,(i-1)*6+5];
    plotPSTH(spktimes, evts, 'kernel', 'gauss', 'eventLims', [x.lick(1) x.lick(end)], ...
        'color', cfg.lickColor, 'subplots', subplots, 'scatSize', scatsz, 'binNum', 100); rasterize
    if i==nunits; xlabel('time from lick (s)'); else; set(gca, 'XTickLabel', []); end
    set(gca, 'xtick', [x.lick(1) 0 x.lick(end)])
    ylabel(''); subplot(nunits*2,3,(i-1)*6+2); ylabel('');
end


% wisk
for i = 1:nunits
    session = units.wisk.session{i};
    unit = units.wisk.unit(i);
    
    e = load(['E:\lab_files\paper2\modelling\neuralData\' session '_neuralData.mat']);
    spktimes = e.spkTimes{e.unit_ids==unit};
    spkrate  = e.spkRates(e.unit_ids==unit,:);
    tstart = e.timeStamps(find(~isnan(spkrate), 1, 'first'));
    tend = e.timeStamps(find(~isnan(spkrate), 1, 'last'));
    
    p = load(['E:\lab_files\paper2\modelling\predictors\' session '_predictors.mat']);
    evts = p.predictors{'whiskerContact', 'data'}{1};
    evts = evts(evts>tstart & evts<tend);
    
%     subplots = [6,ncols,i+4*ncols; 6,ncols,i+5*ncols];
    subplots = [nunits*2,3,(i-1)*6+3; nunits*2,3,(i-1)*6+6];
    plotPSTH(spktimes, evts, 'kernel', 'doubleExp', 'eventLims', [x.wisk(1) x.wisk(end)], ...
        'color', cfg.wiskColor, 'subplots', subplots, 'scatSize', scatsz, 'binNum', 100); rasterize
    set(gca, 'xtick', [x.wisk(1) 0 x.wisk(end)])
    if i==nunits; xlabel('time from whisker contact (s)'); else; set(gca, 'XTickLabel', []); end
    ylabel(''); subplot(nunits*2,3,(i-1)*6+3); ylabel('');
end

saveas(gcf, 'E:\lab_files\paper2\paper_figures\matlab\sample_units', 'svg')


%%  test k stats

close all
for i = 1:10
    smp = unitInfo.step_phase{i};
    [k, p] = ktestnull(smp, 'plot', true, 'res', 40, 'bootstraps', 100);
    pause(.1)
end






































