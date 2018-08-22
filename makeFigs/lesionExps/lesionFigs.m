%% load session info

sessionInfo = readtable([getenv('OBSDATADIR') 'sessions\sessionInfo.xlsx'], 'Sheet', 'lesionNotes');
sessionInfo = sessionInfo(sessionInfo.include==1 & ~cellfun(@isempty, sessionInfo.session),:);



%% get kinematic data
obsPos = -0.0087;
kinData = getKinematicData4(sessionInfo.session, [], obsPos);

% incorporate condition information into kinData struct
for i = 1:length(kinData)
    sessionInfoBin = strcmp(sessionInfo.session, kinData(i).session);
    
    brainRegion = sessionInfo.brainRegion{sessionInfoBin};
    sideOfBrain = sessionInfo.side{sessionInfoBin};
    condition = sessionInfo.preOrPost{sessionInfoBin};
    
    kinData(i).brainRegion = brainRegion;
    kinData(i).sideOfBrain = sideOfBrain;
    kinData(i).condition = condition;
end

data = kinData; save([getenv('OBSDATADIR') 'matlabData\lesionKinematicData.mat'], 'data'); clear data;


%% load kinematic data

load([getenv('OBSDATADIR') 'matlabData\lesionKinematicData.mat'], 'data');
kinData = data; clear data;

%% get speed and avoidance data

speedAvoidanceData = getSpeedAndObsAvoidanceData(sessionInfo.session, false);

% incorporate condition information into kinData struct
for i = 1:length(speedAvoidanceData)
    sessionInfoBin = strcmp(sessionInfo.session, speedAvoidanceData(i).session);
    
    brainRegion = sessionInfo.brainRegion{sessionInfoBin};
    sideOfBrain = sessionInfo.side{sessionInfoBin};
    condition = sessionInfo.preOrPost{sessionInfoBin};
    
    speedAvoidanceData(i).brainRegion = brainRegion;
    speedAvoidanceData(i).sideOfBrain = sideOfBrain;
    speedAvoidanceData(i).condition = condition;
end

data = speedAvoidanceData; save([getenv('OBSDATADIR') 'matlabData\lesionSpeedAvoidanceData.mat'], 'data'); clear data;


%% load speed and avoidance data

load([getenv('OBSDATADIR') 'matlabData\lesionSpeedAvoidanceData.mat'], 'data');
speedAvoidanceData = data; clear data;

%% plot baseline kinematics

mice = unique(sessionInfo.mouse);

for j = 1:length(mice)
    
    mouseBins = strcmp(sessionInfo.mouse, mice{j});
    brainRegion = unique(sessionInfo.brainRegion(mouseBins));
    preBins = strcmp(sessionInfo.preOrPost, 'pre');
    postBins = strcmp(sessionInfo.preOrPost, 'post');
    sides = unique(sessionInfo.side(mouseBins));
    
    for k = 1:length(sides)
        sideBins = strcmp(sessionInfo.side, sides{k});
        preSessions = sessionInfo.session(mouseBins & sideBins & preBins);
        postSessions = sessionInfo.session(mouseBins & sideBins & postBins);
        
        conditionBins = nan(1,length(kinData));
        conditionBins(ismember({kinData.session}, preSessions)) = 1;
        conditionBins(ismember({kinData.session}, postSessions)) = 2;
        
        conditionLabels = {['pre - ' sides{k} ' ' brainRegion{1}], ['post - ' sides{k} ' ' brainRegion{1}]};
        plotBaselineKinematics(kinData, conditionBins, conditionLabels, mice{j})
    end
end


%% plot dv averages

manipulationBarPlots(speedAvoidanceData, 'lesions');
saveas(gcf, [getenv('OBSDATADIR') 'figures\lesions\lesionBarPlots.png']);
savefig([getenv('OBSDATADIR') 'figures\lesions\leionsBarPlots.fig'])


%% plot dvs across sessions

manipulationAcrossSessions(speedAvoidanceData, 'lesions');
saveas(gcf, [getenv('OBSDATADIR') 'figures\lesions\lesionAcrossSessions.png']);
savefig([getenv('OBSDATADIR') 'figures\lesions\lesionAcrossSessions.fig'])


%% paw height by obs height

binNames = {'pre', 'sen lesion', 'mtc lesion'};
% leftSideBins = strcmp({kinData.sideOfBrain}, 'left');
% contraSides = ones(1,length(kinData));
% contraSides(leftSideBins) = 3;
% contraFirstBins = [kinData.firstPawOver] == contraSides;

bins = zeros(1,length(kinData));
bins(strcmp({kinData.condition}, 'pre')) = 1;
bins(strcmp({kinData.condition}, 'post') & strcmp({kinData.brainRegion}, 'sen')) = 2;
bins(strcmp({kinData.condition}, 'post') & strcmp({kinData.brainRegion}, 'mtc') & contraFirstBins) = 3;

scatterObsVsPawHeights(kinData, bins, binNames);
saveas(gcf, fullfile(getenv('OBSDATADIR'), 'figures/lesions/heightShapingScatter.png'));
savefig(fullfile(getenv('OBSDATADIR'), 'figures/lesions/heightShapingScatter.fig'))


%% plot kinematics



% incorporate condition information into kinData struct
for i = 1:length(kinData)
    sessionInfoBin = strcmp(sessionInfo.session, kinData(i).session);
    brainRegion = sessionInfo.brainRegion{sessionInfoBin};
    preOrPost = sessionInfo.preOrPost{sessionInfoBin};
    if strcmp(sessionInfo.side{sessionInfoBin}, 'left')
        contraLimb = 3;
    elseif strcmp(sessionInfo.side{sessionInfoBin}, 'right')
        contraLimb = 2;
    else
        contraLimb = nan;
    end
    
    kinData(i).brainRegion = brainRegion;
    kinData(i).preOrPost = preOrPost;
    kinData(i).contraLimb = contraLimb;
end

% get trial bins

mice = unique(sessionInfo.mouse);

firstPawContraBins = [kinData.firstPawOver]==[kinData.contraLimb];
firstPawIpsiBins = [kinData.firstPawOver]~=[kinData.contraLimb];
preBins = strcmp({kinData.preOrPost}, 'pre');
postBins = strcmp({kinData.preOrPost}, 'post');
lightOnBins = [kinData.isLightOn];
mtcBins = strcmp({kinData.brainRegion}, 'mtc');
senBins = strcmp({kinData.brainRegion}, 'sen');
validBins = ~[kinData.isWheelBreak] & [kinData.isObsAvoided] & ismember({kinData.mouse}, mice);

% % mtc, light on
% bins = zeros(1,length(kinData));
% figBins = mtcBins & validBins;
% bins(firstPawIpsiBins & preBins & figBins) = 1; % ipsi
% bins(firstPawContraBins & preBins & figBins) = 2; % contra
% bins(postBins & figBins) = 3; % vehicle
% binLabels = {'ipsi', 'contra', 'vehicle ipsi'};
% plotObsHeightTrajectories(kinData, bins, binLabels, 'motor cortex lesions, light on')
% 
% % mtc, light off
% bins = zeros(1,length(kinData));
% figBins = mtcBins & ~lightOnBins & validBins;
% bins(firstPawIpsiBins & preBins & figBins) = 1; % ipsi
% bins(firstPawContraBins & preBins & figBins) = 2; % contra
% bins(postBins & figBins) = 3; % vehicle
% binLabels = {'ipsi', 'contra', 'vehicle'};
% plotObsHeightTrajectories(kinData, bins, binLabels, 'motor cortex lesions, light off')

% sen, light on
bins = zeros(1,length(kinData));
figBins = senBins & lightOnBins & validBins;
bins(preBins & figBins) = 1;
bins(postBins & figBins) = 2;
binLabels = {'pre', 'post'};
plotObsHeightTrajectories(kinData, bins, binLabels, 'sensory cortex lesions, light on')

% sen, light off
bins = zeros(1,length(kinData));
figBins = senBins & ~lightOnBins & validBins;
bins(preBins & figBins) = 1;
bins(postBins & figBins) = 2;
binLabels = {'pre', 'post'};
plotObsHeightTrajectories(kinData, bins, binLabels, 'sensory cortex lesions, light off')









