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


%% load kinematic data

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


%% load speed and avoidance data

load([getenv('OBSDATADIR') 'matlabData\sensoryDependenceSpeedAvoidanceData.mat'], 'data');
speedAvoidanceData = data; clear data;

%% plot speed and success

sensoryDependenceSuccessAndSpeed(speedAvoidanceData);
saveas(gcf, fullfile(getenv('OBSDATADIR'), 'figures/sensoryDependence/speedAndSuccess.png'));
savefig(fullfile(getenv('OBSDATADIR'), 'figures/sensoryDependence/speedAndSuccess.fig'))


%% plot sensory dependence kinematics

sensoryDependenceKinematics(kinData);
saveas(gcf, fullfile(getenv('OBSDATADIR'), 'figures/sensoryDependence/sensoryDependenceKinematics.png'));
savefig(fullfile(getenv('OBSDATADIR'), 'figures/sensoryDependence/sensoryDependenceKinematics.fig'))


%% paw height by obs height

wiskSessions = unique({kinData(strcmp({kinData.preOrPost}, 'pre')).session});
noWiskSessions = unique({kinData(strcmp({kinData.preOrPost}, 'post')).session});
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
