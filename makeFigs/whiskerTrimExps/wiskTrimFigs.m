%% load session info

sessionInfo = readtable(fullfile(getenv('OBSDATADIR'), 'spreadSheets', 'experimentMetadata.xlsx'), 'Sheet', 'whiskerTrimNotes');
% sessionInfo = sessionInfo(sessionInfo.include==1 & ~cellfun(@isempty, sessionInfo.session),:);
sessionInfo = sessionInfo(~cellfun(@isempty, sessionInfo.session),:);
sessionInfo.condition = cellfun(@(x,y) [x ' ' y], sessionInfo.whiskerSides, sessionInfo.whiskers, 'UniformOutput', false);


%% compute kinematic data

loadPreviousData = true;

fileName = fullfile(getenv('OBSDATADIR'), 'matlabData', 'whiskerTrimKinematicData.mat');
if loadPreviousData && exist(fileName, 'file')
    load(fileName, 'data');
    kinData = getKinematicData4(sessionInfo.session, sessionInfo, data);
else
    kinData = getKinematicData4(sessionInfo.session, sessionInfo, []);
end
kinData = kinData(~[kinData.isLightOn]); % use only light off trials for these analyses
data = kinData; save(fileName, 'data', '-v7.3', '-nocompression'); clear data;


%% load kinematic data

load(fullfile(getenv('OBSDATADIR'), 'matlabData','whiskerTrimKinematicData.mat'), 'data');
kinData = data; clear data;

%% get speed and avoidance data

speedAvoidanceData = getSpeedAndObsAvoidanceData(sessionInfo.session, sessionInfo, true);
speedAvoidanceData = speedAvoidanceData(~[speedAvoidanceData.isLightOn]); % use only light off trials for these analyses
data = speedAvoidanceData; save(fullfile(getenv('OBSDATADIR'), 'matlabData','whiskerTrimSpeedAvoidanceData.mat'), 'data'); clear data;


%% load speed and avoidance data

load(fullfile(getenv('OBSDATADIR'), 'matlabData','whiskerTrimSpeedAvoidanceData.mat'), 'data');
speedAvoidanceData = data; clear data;


%%  compute dependent measures

dvs = {'success', 'forePawErrRate', 'speed', 'bodyAngleContra', 'contraFirstRate', 'bigStepProb', 'pawHgt', 'hgtShaping'};
sessionDvs = getSessionDvs(dvs, speedAvoidanceData, kinData);
disp('finished computing dependent measures!')
% sessionDvs = getSessionDvs(dvs, speedAvoidanceData);

%% bar plots

dvs = {'success', 'speed', 'bodyAngleContra', 'forePawErrRateIpsi', 'forePawErrRateContra', ...
    'contraFirstRate', 'bigStepProbIpsi', 'bigStepProbContra', 'pawHgtIpsi', 'pawHgtContra', ...
    'hgtShapingIpsi', 'hgtShapingContra'};
barPlots(sessionDvs, dvs, 'whiskerTrimming')
saveas(gcf, fullfile(getenv('OBSDATADIR'), 'figures/', 'whiskerTrim', '/whiskerTrimBarPlots.png'));
savefig(fullfile(getenv('OBSDATADIR'), 'figures/', 'whiskerTrim', '/whiskerTrimBarPlots.fig'))
%% sessions over time plots

% miceToShow = {'den17', 'den18', 'den19'}; % set to 'all' to show all mice
miceToShow = {'den10', 'den12'}; % set to 'all' to show all mice

if strcmp(miceToShow, 'all'); bins = true(1,length(sessionDvs)); else; bins = ismember({sessionDvs.mouse}, miceToShow); end
plotAcrossSessions(sessionDvs(bins), dvs, 'whiskerTrimming')
saveas(gcf, fullfile(getenv('OBSDATADIR'), 'figures/', 'whiskerTrim', '/whiskerTrimAcrossSessions.png'));
savefig(fullfile(getenv('OBSDATADIR'), 'figures/', 'whiskerTrim', '/whiskerTrimAcrossSessions.fig'))
%% paw height by obs height

binNames = {'ipsi', 'contra'};
conditions = {'bilateral full', 'unilateral full'};
conditionNames = {'bilateral full (ipsi)', 'bilateral full (contra)', 'unilateral full (ipsi)', 'unilateral full (contra)'};

conditionBins = zeros(1,length(kinData));
conditionBins(strcmp({kinData.condition}, conditions{1}) & [kinData.ipsiPawFirst]) = 1;
conditionBins(strcmp({kinData.condition}, conditions{1}) & [kinData.contraPawFirst]) = 2;
conditionBins(strcmp({kinData.condition}, conditions{2}) & [kinData.ipsiPawFirst]) = 3;
conditionBins(strcmp({kinData.condition}, conditions{2}) & [kinData.contraPawFirst]) = 4;

% conditionBins([kinData.conditionNum]<3 & ~strcmp({kinData.mouse}, 'den10')) = 0; % restrict trials
conditionBins([kinData.conditionNum]==1) = 0; % restrict trials
% conditionBins(~strcmp({kinData.mouse}, 'den10')) = 0; % restrict trialsr

scatterObsVsPawHeights(kinData, conditionBins, conditionNames);
saveas(gcf, fullfile(getenv('OBSDATADIR'), 'figures/', 'whiskerTrim', '/heightShapingScatter.png'));
savefig(fullfile(getenv('OBSDATADIR'), 'figures/', 'whiskerTrim', '/heightShapingScatter.fig'))

