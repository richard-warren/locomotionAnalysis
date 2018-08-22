%% load session info

sessionInfo = readtable([getenv('OBSDATADIR') 'sessions\sessionInfo.xlsx'], 'Sheet', 'whiskerTrimNotes');
sessionInfo = sessionInfo(sessionInfo.include==1 & ~cellfun(@isempty, sessionInfo.session),:);



%% get kinematic data
obsPos = -0.0087;
kinData = getKinematicData4(sessionInfo.session, [], obsPos);

% incorporate condition information into kinData struct
for i = 1:length(kinData)
    sessionInfoBin = strcmp(sessionInfo.session, kinData(i).session);
    preOrPost = sessionInfo.preOrPost{sessionInfoBin};
    kinData(i).preOrPost = preOrPost;
end

data = kinData; save([getenv('OBSDATADIR') 'matlabData\sensoryDependenceKinematicData.mat'], 'data'); clear data;


%% load previous kinematic data

load([getenv('OBSDATADIR') 'matlabData\sensoryDependenceKinematicData.mat'], 'data');
kinData = data; clear data;

%% get speed and avoidance data

speedAvoidanceData = getSpeedAndObsAvoidanceData(sessionInfo.session, false);

% incorporate condition information into kinData struct
for i = 1:length(speedAvoidanceData)
    sessionInfoBin = strcmp(sessionInfo.session, speedAvoidanceData(i).session);
    preOrPost = sessionInfo.preOrPost{sessionInfoBin};
    speedAvoidanceData(i).preOrPost = preOrPost;
end

data = speedAvoidanceData; save([getenv('OBSDATADIR') 'matlabData\sensoryDependenceSpeedAvoidanceData.mat'], 'data'); clear data;


%% plot speed and success

sensoryDependenceSuccessAndSpeed(speedAvoidanceData);
saveas(gcf, fullfile(getenv('OBSDATADIR'), 'figures/sensoryDependence/speedAndSuccess.png'));
savefig(fullfile(getenv('OBSDATADIR'), 'figures/sensoryDependence/speedAndSuccess.fig'))


%% plot sensory dependence kinematics

sensoryDependenceKinematics(kinData);
saveas(gcf, fullfile(getenv('OBSDATADIR'), 'figures/sensoryDependence/sensoryDependenceKinematics.png'));
savefig(fullfile(getenv('OBSDATADIR'), 'figures/sensoryDependence/sensoryDependenceKinematics.fig'))


%% paw height by obs height

wiskSessions = unique(sessionInfo.session(strcmp(sessionInfo.preOrPost, 'pre')));
noWiskSessions = unique(sessionInfo.session(strcmp(sessionInfo.preOrPost, 'post')));
wiskLightBins = ismember({kinData.session}, wiskSessions) & [kinData.isLightOn];
wiskBins = ismember({kinData.session}, wiskSessions) & ~[kinData.isLightOn];
lightBins = ismember({kinData.session}, noWiskSessions) & [kinData.isLightOn];
neitherBins = ismember({kinData.session}, noWiskSessions) & ~[kinData.isLightOn];
conditionBins = {wiskLightBins, wiskBins, lightBins, neitherBins};
binNames = {'wisk+light', 'wisk', 'light', 'neither'};
bins = zeros(1,length(kinData));
for i = 1:4; bins(conditionBins{i}) = i; end

scatterObsVsPawHeights(kinData, bins, binNames);
saveas(gcf, fullfile(getenv('OBSDATADIR'), 'figures/sensoryDependence/heightShapingScatter.png'));
savefig(fullfile(getenv('OBSDATADIR'), 'figures/sensoryDependence/heightShapingScatter.fig'))


%% plot baseline kinematics

mice = unique(sessionInfo.mouse);

for j = 1:length(mice)
    
    mouseBins = strcmp(sessionInfo.mouse, mice{j});
    preBins = strcmp(sessionInfo.preOrPost, 'pre');
    postBins = strcmp(sessionInfo.preOrPost, 'post');
    
    preSessions = sessionInfo.session(mouseBins & preBins);
    postSessions = sessionInfo.session(mouseBins & postBins);

    conditionBins = nan(1,length(kinData));
    conditionBins(ismember({kinData.session}, preSessions)) = 1;
    conditionBins(ismember({kinData.session}, postSessions)) = 2;

    conditionLabels = {'pre trimming', 'post trimming'};
    plotBaselineKinematics(kinData, conditionBins, conditionLabels, mice{j})
end


%% plot success and speed averages


% saveas(gcf, [getenv('OBSDATADIR') 'figures\sensoryDependencePlots.png']);
% savefig([getenv('OBSDATADIR') 'figures\sensoryDependenceBarPlots.fig'])


