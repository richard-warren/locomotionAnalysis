%% load session info

sessionInfo = readtable([getenv('OBSDATADIR') 'sessions\sessionInfo.xlsx'], 'Sheet', 'muscimolNotes');
sessionInfo = sessionInfo(sessionInfo.include==1 & ~cellfun(@isempty, sessionInfo.session),:);



%% get kinematic data
kinData = getKinematicData4(sessionInfo.session, []);

% incorporate condition information into kinData struct
for i = 1:length(kinData)
    sessionInfoBin = strcmp(sessionInfo.session, kinData(i).session);
    
    brainRegion = sessionInfo.brainRegion{sessionInfoBin};
    sideOfBrain = sessionInfo.side{sessionInfoBin};
    condition = sessionInfo.injectedSubstance{sessionInfoBin};
    
    kinData(i).brainRegion = brainRegion;
    kinData(i).sideOfBrain = sideOfBrain;
    kinData(i).condition = condition;
end

data = kinData; save([getenv('OBSDATADIR') 'matlabData\muscimolKinematicData.mat'], 'data'); clear data;


%% load kinematic data

load([getenv('OBSDATADIR') 'matlabData\muscimolKinematicData.mat'], 'data');
kinData = data; clear data;
disp('muscimol kinematic data loaded!')

%% get speed and avoidance data

speedAvoidanceData = getSpeedAndObsAvoidanceData(sessionInfo.session, false);

% incorporate condition information into kinData struct
for i = 1:length(speedAvoidanceData)
    sessionInfoBin = strcmp(sessionInfo.session, speedAvoidanceData(i).session);
    
    brainRegion = sessionInfo.brainRegion{sessionInfoBin};
    sideOfBrain = sessionInfo.side{sessionInfoBin};
    condition = sessionInfo.injectedSubstance{sessionInfoBin};
    
    speedAvoidanceData(i).brainRegion = brainRegion;
    speedAvoidanceData(i).sideOfBrain = sideOfBrain;
    speedAvoidanceData(i).condition = condition;
end

data = speedAvoidanceData; save([getenv('OBSDATADIR') 'matlabData\muscimolSpeedAvoidanceData.mat'], 'data'); clear data;


%% load speed and avoidance data

load([getenv('OBSDATADIR') 'matlabData\muscimolSpeedAvoidanceData.mat'], 'data');
speedAvoidanceData = data; clear data;
disp('muscimol speed avoidance data loaded!')

%% plot baseline kinematics

mice = unique(sessionInfo.mouse);

for j = 1:length(mice)
    
    mouseBins = strcmp(sessionInfo.mouse, mice{j});
    brainRegion = unique(sessionInfo.brainRegion(mouseBins));
    musBins = strcmp(sessionInfo.injectedSubstance, 'muscimol');
    vehBins = strcmp(sessionInfo.injectedSubstance, 'saline');
    sides = unique(sessionInfo.side(mouseBins & (musBins | vehBins)));
    
    for k = 1:length(sides)
        sideBins = strcmp(sessionInfo.side, sides{k});
        musSessions = sessionInfo.session(mouseBins & sideBins & musBins);
        vehSessions = sessionInfo.session(mouseBins & sideBins & vehBins);
        
        conditionBins = nan(1,length(kinData));
        conditionBins(ismember({kinData.session}, musSessions)) = 1;
        conditionBins(ismember({kinData.session}, vehSessions)) = 2;
        
        conditionLabels = {['muscimol - ' sides{k} ' ' brainRegion{1}], ['vehicle - ' sides{k} ' ' brainRegion{1}]};
        plotBaselineKinematics(kinData, conditionBins, conditionLabels, mice{j})
    end
end


%% plot dv averages

manipulationBarPlots(speedAvoidanceData, {'saline', 'muscimol'}, 'muscimol');
saveas(gcf, [getenv('OBSDATADIR') 'figures\muscimol\muscimolBarPlots.png']);
savefig([getenv('OBSDATADIR') 'figures\muscimol\muscimolBarPlots.fig'])

%% plot dvs across sessions

manipulationAcrossSessions(speedAvoidanceData, {'saline', 'muscimol'}, 'muscimol');
saveas(gcf, [getenv('OBSDATADIR') 'figures\muscimol\muscimolAcrossSessions.png']);
savefig([getenv('OBSDATADIR') 'figures\muscimol\muscimolAcrossSessions.fig'])


%% paw height by obs height

binNames = {'saline', 'sen muscimol', 'mtc muscimol ipsi', 'mtc muscimol contra'};
leftSideBins = strcmp({kinData.sideOfBrain}, 'left');
contraSides = ones(1,length(kinData))*2;
contraSides(leftSideBins) = 3;
contraFirstBins = [kinData.firstPawOver] == contraSides;

bins = zeros(1,length(kinData));
bins(strcmp({kinData.condition}, 'saline')) = 1;
bins(strcmp({kinData.condition}, 'muscimol') & strcmp({kinData.brainRegion}, 'sen')) = 2;
bins(strcmp({kinData.condition}, 'muscimol') & strcmp({kinData.brainRegion}, 'mtc') & ~contraFirstBins) = 3;
bins(strcmp({kinData.condition}, 'muscimol') & strcmp({kinData.brainRegion}, 'mtc') & contraFirstBins) = 4;

scatterObsVsPawHeights(kinData, bins, binNames);
saveas(gcf, fullfile(getenv('OBSDATADIR'), 'figures/muscimol/heightShapingScatter.png'));
savefig(fullfile(getenv('OBSDATADIR'), 'figures/muscimol/heightShapingScatter.fig'))


%% plot dvs over time within sessions

% settings
trialBinSize = 40;
conditions = {'saline', 'muscimol'};
brainRegions = {'sen', 'mtc'};
dvs = {'success rate', ...
       'speed (m/s)', ...
       {'body angle towards contra', '(or right) side'}};
dvYLims = [0 1; 0 .8; -15 15];
xMax = 240;


% initializations
figure('Color', 'white', 'MenuBar', 'none', 'Position', [-1396 200 800 600])
dims = [length(dvs), length(brainRegions)]; % subplot dimensions
colors = winter(length(conditions));

% loop over brain regions
for i = 1:length(brainRegions)
    
    % get max trial num
    brainRegionBins = ismember({speedAvoidanceData.session}, sessions);
    maxTrial = max([speedAvoidanceData(brainRegionBins).trialNum]);
    maxTrial = min(maxTrial, xMax);
    maxTrial = floor(maxTrial/trialBinSize) * trialBinSize;
    trialBinEdges = 0:trialBinSize:maxTrial;
    xInds = trialBinEdges(2:end); % the x positions will be centered within the bins
    
    
    for j = 1:length(conditions)
        
        sessions = unique(sessionInfo.session(strcmp(sessionInfo.brainRegion, brainRegions{i}) ...
                                            & strcmp(sessionInfo.injectedSubstance, conditions{j})));
        
        % containers for bins averages across all sessions
        successes = nan(length(sessions), length(trialBinEdges)-1);
        speeds = nan(length(sessions), length(trialBinEdges)-1);
        contraBodyAngles = nan(length(sessions), length(trialBinEdges)-1);
        
        for k = 1:length(sessions)
            
            load([getenv('OBSDATADIR') 'sessions\' sessions{k} '\runAnalyzed.mat'], ...
                'bodyAngles', 'obsOnTimes', 'obsOffTimes', 'frameTimeStamps');
            load([getenv('OBSDATADIR') 'sessions\' sessions{k} '\run.mat'], 'breaks');
            bodyAngles = getTrialBodyAngles(bodyAngles, obsOnTimes, obsOffTimes, frameTimeStamps, breaks);
            sideOfBrain = sessionInfo.side{strcmp(sessionInfo.session, sessions{k})};
            if strcmp(sideOfBrain, 'left'); bodyAngles = -bodyAngles; end
            
            
            for m = 1:length(trialBinEdges)-1
            
                trialBinBins = strcmp({speedAvoidanceData.session}, sessions{k}) & ...
                               [speedAvoidanceData.trialNum]>trialBinEdges(m) & ...
                               [speedAvoidanceData.trialNum]<trialBinEdges(m+1);

                % get speed, success rate, and body angle
                successes(k,m) = mean([speedAvoidanceData(trialBinBins).isObsAvoided]);
                speeds(k,m) = mean([speedAvoidanceData(trialBinBins).avgVel]);
                bodyAngleBins = 1:length(bodyAngles) > trialBinEdges(m) & 1:length(bodyAngles) < trialBinEdges(m+1);
                contraBodyAngles(k,m) = nanmean(bodyAngles(bodyAngleBins));
            end
            
            
            % plot session data
            % plot session success
            subplot(dims(1), dims(2), i)
            line(xInds, successes(k,:), 'color', colors(j,:)); hold on

            % plot session speed
            subplot(dims(1), dims(2), i+1*length(brainRegions))
            line(xInds, speeds(k,:), 'color', colors(j,:)); hold on

            % plot session body angle
            subplot(dims(1), dims(2), i+2*length(brainRegions))
            line(xInds, contraBodyAngles(k,:), 'color', colors(j,:)); hold on
        end
        
        
        % plot condition means
        
        % mean successes
        subplot(dims(1), dims(2), i)
        line(xInds, nanmean(successes,1), 'color', colors(j,:), 'linewidth', 4); hold on

        % mean speeds
        subplot(dims(1), dims(2), i+1*length(brainRegions))
        line(xInds, nanmean(speeds,1), 'color', colors(j,:), 'linewidth', 4); hold on
        
        % mean session body angle
        subplot(dims(1), dims(2), i+2*length(brainRegions))
        line(xInds, nanmean(contraBodyAngles,1), 'color', colors(j,:), 'linewidth', 4); hold on
    end
    
    % add condition labels
    xs = linspace(trialBinEdges(2)*1.2, trialBinEdges(end)*.8, length(conditions));
    for j = 1:length(conditions)
        text(xs(j), dvYLims(end,1)+(dvYLims(end,2)-dvYLims(end,1))*.2, conditions{j}, 'Color', colors(j,:));
    end
end
    


% pimp figs
ind = 1;
for i = 1:length(dvs)
    for j = 1:length(brainRegions)
        subplot(dims(1), dims(2), ind);
        
        set(gca, 'xlim', [trialBinEdges(2) trialBinEdges(end)], 'xtick', trialBinEdges(2:end), 'YLim', dvYLims(i,:));
        if i==1; title(brainRegions{j}); end
        if j==1; ylabel(dvs{i}); end
        ind = ind+1;
    end
end


saveas(gcf, [getenv('OBSDATADIR') 'figures\muscimolSessionProgressPlots.png']);
savefig([getenv('OBSDATADIR') 'figures\muscimolSessionProgressPlotsPlots.fig'])




%% plot kinematics



% incorporate condition information into kinData struct
for i = 1:length(kinData)
    sessionInfoBin = strcmp(sessionInfo.session, kinData(i).session);
    brainRegion = sessionInfo.brainRegion{sessionInfoBin};
    injectedSubstance = sessionInfo.injectedSubstance{sessionInfoBin};
    if strcmp(sessionInfo.side{sessionInfoBin}, 'left')
        contraLimb = 3;
    elseif strcmp(sessionInfo.side{sessionInfoBin}, 'right')
        contraLimb = 2;
    else
        contraLimb = nan;
    end
    
    kinData(i).brainRegion = brainRegion;
    kinData(i).injectedSubstance = injectedSubstance;
    kinData(i).contraLimb = contraLimb;
end

% get trial bins

mice = unique(sessionInfo.mouse);

firstPawContraBins = [kinData.firstPawOver]==[kinData.contraLimb];
firstPawIpsiBins = [kinData.firstPawOver]~=[kinData.contraLimb];
musBins = strcmp({kinData.injectedSubstance}, 'muscimol');
vehBins = strcmp({kinData.injectedSubstance}, 'saline');
lightOnBins = [kinData.isLightOn];
mtcBins = strcmp({kinData.brainRegion}, 'mtc');
senBins = strcmp({kinData.brainRegion}, 'sen');
% validBins = ~[kinData.isWheelBreak] & [kinData.isObsAvoided] & ismember({kinData.mouse}, mice);
validBins = ~[kinData.isWheelBreak] & ~[kinData.isObsAvoided] & ismember({kinData.mouse}, mice);

% mtc, light on
bins = zeros(1,length(kinData));
figBins = mtcBins & validBins;
bins(firstPawIpsiBins & musBins & figBins) = 1; % ipsi
bins(firstPawContraBins & musBins & figBins) = 2; % contra
bins(vehBins & figBins) = 3; % vehicle
binLabels = {'ipsi', 'contra', 'vehicle ipsi'};
plotObsHeightTrajectories(kinData, bins, binLabels, 'motor cortex lesions, light on')

% mtc, light off
bins = zeros(1,length(kinData));
figBins = mtcBins & ~lightOnBins & validBins;
bins(firstPawIpsiBins & musBins & figBins) = 1; % ipsi
bins(firstPawContraBins & musBins & figBins) = 2; % contra
bins(vehBins & figBins) = 3; % vehicle
binLabels = {'ipsi', 'contra', 'vehicle'};
plotObsHeightTrajectories(kinData, bins, binLabels, 'motor cortex lesions, light off')

% % sen, light on
% bins = zeros(1,length(kinData));
% figBins = senBins & lightOnBins & validBins;
% bins(musBins & figBins) = 1; % mus
% bins(vehBins & figBins) = 2; % vehicle
% binLabels = {'muscimol', 'vehicle'};
% plotObsHeightTrajectories(kinData, bins, binLabels, 'sensory cortex lesions, light on')
% 
% % sen, light off
% bins = zeros(1,length(kinData));
% figBins = senBins & ~lightOnBins & validBins;
% bins(musBins & figBins) = 1; % mus
% bins(vehBins & figBins) = 2; % vehicle
% binLabels = {'muscimol', 'vehicle'};
% plotObsHeightTrajectories(kinData, bins, binLabels, 'sensory cortex lesions, light off')















