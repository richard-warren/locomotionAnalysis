%% load session info

sessionInfo = readtable([getenv('OBSDATADIR') 'sessions\sessionInfo.xlsx'], 'Sheet', 'whiskerTrimNotes');
sessionInfo = sessionInfo(sessionInfo.include==1 & ~cellfun(@isempty, sessionInfo.session),:);


%% get kinematic data

wiskBins = strcmp(sessionInfo.preOrPost, 'pre');
noWiskBins = strcmp(sessionInfo.preOrPost, 'post');

kinDataWisk = getKinematicData4(sessionInfo.session(wiskBins), sessionInfo, []);
obsPos = nanmedian([kinDataWisk.obsPos]);
kinDataNoWisk = getKinematicData4(sessionInfo.session(noWiskBins), sessionInfo, [], obsPos);
kinData = cat(2, kinDataWisk, kinDataNoWisk);

data = kinData; save([getenv('OBSDATADIR') 'matlabData\sensoryDependenceKinematicData.mat'], 'data');
clear data kinDataWisk kinDataNoWisk;

%% load kinematic data

load([getenv('OBSDATADIR') 'matlabData\sensoryDependenceKinematicData.mat'], 'data');
kinData = data; clear data;

%% get speed and avoidance data

speedAvoidanceData = getSpeedAndObsAvoidanceData(sessionInfo.session, sessionInfo, false);
data = speedAvoidanceData; save([getenv('OBSDATADIR') 'matlabData\sensoryDependenceSpeedAvoidanceData.mat'], 'data'); clear data;


%% load speed and avoidance data

load([getenv('OBSDATADIR') 'matlabData\sensoryDependenceSpeedAvoidanceData.mat'], 'data');
speedAvoidanceData = data; clear data;

%% plot speed and success

sensoryDependenceSuccessAndSpeed(speedAvoidanceData);
saveas(gcf, fullfile(getenv('OBSDATADIR'), 'figures/sensoryDependence/speedAndSuccess.png'));
savefig(fullfile(getenv('OBSDATADIR'), 'figures/sensoryDependence/speedAndSuccess.fig'))


%% plot sensory dependence kinematics

bins = zeros(1,length(kinData));
binNames = {'wisk+light', 'wisk', 'light', 'neither'};
bins(strcmp({kinData.preOrPost}, 'pre') & [kinData.isLightOn]) = 1;
bins(strcmp({kinData.preOrPost}, 'pre') & ~[kinData.isLightOn]) = 2;
bins(strcmp({kinData.preOrPost}, 'post') & [kinData.isLightOn]) = 3;
bins(strcmp({kinData.preOrPost}, 'post') & ~[kinData.isLightOn]) = 4;


close all; plotObsHeightTrajectories(kinData, bins, binNames, 'sensoryDependence')

% sensoryDependenceKinematics(kinData);
% saveas(gcf, fullfile(getenv('OBSDATADIR'), 'figures/sensoryDependence/sensoryDependenceKinematics.png'));
% savefig(fullfile(getenv('OBSDATADIR'), 'figures/sensoryDependence/sensoryDependenceKinematics.fig'))


%% paw height by obs height

bins = zeros(1,length(kinData));
binNames = {'wisk+light', 'wisk', 'light', 'neither'};
bins(strcmp({kinData.preOrPost}, 'pre') & [kinData.isLightOn]) = 1;
bins(strcmp({kinData.preOrPost}, 'pre') & ~[kinData.isLightOn]) = 2;
bins(strcmp({kinData.preOrPost}, 'post') & [kinData.isLightOn]) = 3;
bins(strcmp({kinData.preOrPost}, 'post') & ~[kinData.isLightOn]) = 4;

scatterObsVsPawHeights(kinData, bins, binNames);
saveas(gcf, fullfile(getenv('OBSDATADIR'), 'figures/sensoryDependence/heightShapingScatter.png'));
savefig(fullfile(getenv('OBSDATADIR'), 'figures/sensoryDependence/heightShapingScatter.fig'))


%% make sensory dependence vids

dir = 'C:\Users\rick\Desktop\';
titles = {'wisk and light', 'wisk', 'light', 'neither'};
sessions = {'180731_007', '180801_005'}; % first session pre, second session post wisk trimming
trialsToShow = 15;
titleInd = 1;

for session = 1:2
    load([getenv('OBSDATADIR') 'sessions\' sessions{session} '\runAnalyzed.mat'], 'isLightOn');
    for condition = 1:2
        if condition==1; trials = find(isLightOn); else; trials = find(~isLightOn); end
        trials = trials([1:trialsToShow] + round(length(trials)/2)); % take the trials from the middle of the session
        makeVidWisk(fullfile(dir, titles{titleInd}), sessions{session}, [-.05 .1], .1, trials);
        titleInd = titleInd+1;
    end
end


