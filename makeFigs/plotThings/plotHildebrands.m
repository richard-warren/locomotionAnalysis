function plotHildebrands(sessionInfo, varargin)



% settings
s.refPaw = 2;  % 2 or 3 // hildebrand plots start with refPaw swing and end with subsequent swing
s.bins = 200; % length of x grid for hildebrand plots
s.colors = hsv(4);
s.plotMice = false;  % whether to plot all mice in addition to group means
s.stepPercentiles = [5 95];  % for overal plots, only include steps with durations within these percetile limits

s.pawNames = {'left hind', 'left fore', 'right fore', 'right hind'};
s.stepTypes = {'leading hind', 'leading fore', 'trailing fore', 'trailing hind'};  % assuming s.refPaw=2


% initializations
if exist('varargin', 'var'); for i = 1:2:length(varargin); s.(varargin{i}) = varargin{i+1}; end; end  % reassign settings passed in varargin
c = permute(repmat(s.colors,1,1,s.bins), [1,3,2]);  % matrix of step type colors
mice = unique(sessionInfo.mouse);
if s.refPaw~=2; s.stepTypes = s.stepTypes([4 3 2 1]); end




[H_all, H_obs] = deal(nan(length(mice), 4, s.bins));  % (mouse) X (paw) X (bins)

for i = 1:length(mice)
    sessions = sessionInfo.session(strcmp(sessionInfo.mouse, mice{i}));
    
    [H_allSub, H_obsSub] = deal(nan(length(sessions), 4, s.bins));  % (mouse) X (paw) X (bins)
    
    for j = 1:length(sessions)
        [H_allSub(j,:,:), H_obsSub(j,:,:)] = hildebrand(sessions{j});
    end
    H_all(i,:,:) = nanmean(H_allSub,1);
    H_obs(i,:,:) = nanmean(H_obsSub,1);
end
 

if s.plotMice; hgt = 100*(length(mice)+1); else; hgt = 100; end
cols = 1 + length(mice)*s.plotMice;
figure('Color', 'white', 'MenuBar', 'none', 'Position', [1975.00 10 632.00 hgt]);



% plot means

% overall
subplot(cols, 2, 1)
colormap gray
imagesc(squeeze(nanmean(H_all,1)))
set(gca, 'Visible', 'on', 'XTick', [], 'YTick', 1:4, 'YTickLabel', s.pawNames, 'TickLength', [0 0])
title('overall')
xlabel('fraction stride')

% hurdling
subplot(cols, 2, 2)
image(squeeze(nanmean(H_obs,1)).*c)
set(gca, 'Visible', 'on', 'XTick', [], 'YTick', 1:4, 'YTickLabel', s.stepTypes, 'TickLength', [0 0])
title('hurdle clearance')
xlabel('fraction stride')


% plot individual mice
if s.plotMice
    for i = 1:length(mice)
        
        % overall
        subplot(cols, 2, 2+(i-1)*2+1)
        colormap gray
        imagesc(squeeze(H_all(i,:,:)))
        set(gca, 'Visible', 'on', 'YTick', [], 'XTick', [])
        ylabel(mice{i})
        
        % hurdling
        subplot(cols, 2, 2+(i-1)*2+2)
        image(squeeze(H_obs(i,:,:)).*c)
        set(gca, 'Visible', 'on', 'YTick', [], 'XTick', [])
    end
end



    
function [hAll, hObs] = hildebrand(session)

    % get session hildebrands for entire session (hAll) and for steps over the hurdle (hObs)
    
    load(fullfile(getenv('OBSDATADIR'), 'sessions', session, 'kinData.mat'), ...
        'kinData', 'stanceBins');
    fprintf('%s: analyzing session...\n', session)

    % overall hildebrand
    isSwing = ~stanceBins;
    swingStarts = find(diff(isSwing(:,s.refPaw))==1)+1;
    
    % limit to swings with valid durations
    swingDurations = diff(swingStarts);
    durationLims = prctile(swingDurations, s.stepPercentiles);
    swingStartInds = find(swingDurations>=durationLims(1) & swingDurations<=durationLims(2))';
    
    hAll = nan(length(swingStarts), 4, s.bins);  % (step number) X (paw) X (position within swing)
    for k = swingStartInds
        inds = swingStarts(k):swingStarts(k+1)-1;
        hAll(k,:,:) = interp1(1:length(inds), double(isSwing(inds,:)), linspace(1,length(inds),s.bins))';
    end
    hAll = squeeze(nanmean(1-hAll,1)); % average and switch to probability of being in stance

    
    % step over hurdle hildebrand
    hObs = nan(length(kinData), 4, s.bins);  % (step number) X (paw) X (position within swing)

    for k = find([kinData.isTrialAnalyzed])
        p = kinData(k).pawOverSequence(1);  % first paw to step over hurdle
        isSwing = ~kinData(k).stanceBins;

        swingStarts = find(diff(isSwing(:,p))==1)+1;
        ind = find(diff(~isnan(kinData(k).modifiedStepIdentities(:,p)))==1, 1, 'last')+1;  % index of start of step over hurdle
        inds = ind : (swingStarts(find(swingStarts==ind)+1)-1);  % inds within trial for step

        hObs(k,:,:) = interp1(1:length(inds), double(isSwing(inds,:)), linspace(1,length(inds),s.bins))';
        if p ~= s.refPaw; hObs(k,:,:) = hObs(k,[4 3 2 1],:); end  % flip left and right sides of body if necessary
    end
    
    hObs = squeeze(nanmean(1-hObs,1)); % average and switch to probability of being in stance
end
end




