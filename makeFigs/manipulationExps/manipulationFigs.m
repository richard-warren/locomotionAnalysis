%% load session info


% settings
manipulation = 'lesion';
maxLesionSession = 4;


if strcmp(manipulation, 'muscimol'); conditions = {'saline', 'muscimol'};
elseif strcmp(manipulation, 'lesion'); conditions = {'pre', 'post'}; end
% sessionInfo = readtable([getenv('OBSDATADIR') 'sessions\sessionInfo.xlsx'], 'Sheet', [manipulation 'Notes']);
sessionInfo = sessionInfo(~cellfun(@isempty, sessionInfo.session) & strcmp(sessionInfo.brainRegion, 'mtc'),:);
sessionInfo = sessionInfo(~cellfun(@isempty, sessionInfo.session),:);



%% compute kinematic data

loadPreviousData = true;

if loadPreviousData
    load([getenv('OBSDATADIR') 'matlabData\' manipulation 'kinematicData.mat'], 'data');
    kinData = getKinematicData4(sessionInfo.session, sessionInfo, data);
else
    kinData = getKinematicData4(sessionInfo.session, sessionInfo, []);
end
data = kinData; save([getenv('OBSDATADIR') 'matlabData\' manipulation 'KinematicData.mat'], 'data', '-v7.3', '-nocompression'); clear data;



%% load kinematic data

load([getenv('OBSDATADIR') 'matlabData\' manipulation 'KinematicData.mat'], 'data');
kinData = data; clear data;
disp([manipulation ' kinematic data loaded!'])


%% compute speed and avoidance data

speedAvoidanceData = getSpeedAndObsAvoidanceData(sessionInfo.session, sessionInfo, false);
data = speedAvoidanceData; save([getenv('OBSDATADIR') 'matlabData\' manipulation 'SpeedAvoidanceData.mat'], 'data'); clear data;


%% load speed and avoidance data

load([getenv('OBSDATADIR') 'matlabData\' manipulation 'SpeedAvoidanceData.mat'], 'data');
speedAvoidanceData = data; clear data;
disp([manipulation ' speed avoidance data loaded!'])

%% ----------

% PLOT THINGS

%  ----------

%% plot dv averages

% don't include too many post lesion sessions if manipulation=='lesion'
if strcmp(manipulation, 'lesion'); includeTrial = [speedAvoidanceData.conditionNum]<=maxLesionSession;
else; includeTrial = true(1,length(speedAvoidanceData)); end

manipulationBarPlots(speedAvoidanceData(includeTrial), conditions, manipulation);
saveas(gcf, [getenv('OBSDATADIR') 'figures\' manipulation '\' manipulation 'BarPlots.png']);
savefig([getenv('OBSDATADIR') 'figures\' manipulation '\' manipulation 'BarPlots.fig'])


%% plot dvs across sessions

manipulationAcrossSessions(speedAvoidanceData, conditions, manipulation);
saveas(gcf, [getenv('OBSDATADIR') 'figures\' manipulation '\' manipulation 'AcrossSession.png']);
savefig([getenv('OBSDATADIR') 'figures\' manipulation '\' manipulation 'AcrossSessions.fig'])

%% paw height by obs height for mtc

binNames = {'ipsi', 'contra', conditions{1}};

bins = zeros(1,length(kinData));
bins(strcmp({kinData.condition}, conditions{2}) & strcmp({kinData.brainRegion}, 'mtc') & [kinData.ipsiPawFirst]) = 1;
bins(strcmp({kinData.condition}, conditions{2}) & strcmp({kinData.brainRegion}, 'mtc') & [kinData.contraPawFirst]) = 2;
bins(strcmp({kinData.condition}, conditions{1})) = 3;

% don't include too many post lesion sessions if manipulation=='lesion'
if strcmp(manipulation, 'lesion'); includeTrial = [kinData.conditionNum]<=maxLesionSession;
else; includeTrial = ones(1,length(kinData)); end

scatterObsVsPawHeights(kinData, bins.*includeTrial, binNames);
saveas(gcf, fullfile(getenv('OBSDATADIR'), 'figures/', manipulation, '/heightShapingScatterMtc.png'));
savefig(fullfile(getenv('OBSDATADIR'), 'figures/', manipulation, '/heightShapingScatterMtc.fig'))


%% plot kinematics

% settings
miceToExclude = {'sen5'};

% get trial bins
contraFirstBins = [kinData.contraPawFirst];
ipsiFirstBins = [kinData.ipsiPawFirst];
controlBins = strcmp({kinData.condition}, conditions{1});
manipBins = strcmp({kinData.condition}, conditions{2});
lightOnBins = [kinData.isLightOn];
mtcBins = strcmp({kinData.brainRegion}, 'mtc');
senBins = strcmp({kinData.brainRegion}, 'sen');

% don't include too many post lesion sessions if manipulation=='lesion'
if strcmp(manipulation, 'lesion'); includeTrial = [kinData.conditionNum]<=maxLesionSession;
else; includeTrial = ones(1,length(kinData)); end

% mtc
bins = zeros(1,length(kinData));
bins(mtcBins & ipsiFirstBins & manipBins) = 1; % ipsi
bins(mtcBins & contraFirstBins & manipBins) = 2; % contra
bins(mtcBins & controlBins) = 3; % pre
binLabels = {'ipsi', 'contra', 'pre'};
plotObsHeightTrajectories(kinData, bins.*includeTrial, binLabels, ['motor cortex ' manipulation])
saveas(gcf, fullfile(getenv('OBSDATADIR'), 'figures', manipulation, 'MtcKinematics.png'));
savefig(fullfile(getenv('OBSDATADIR'), 'figures', manipulation, 'MtcKinematics.fig'))

% sen
bins = zeros(1,length(kinData));
bins(senBins & controlBins) = 1; % control
bins(senBins & manipBins) = 2; % manipluation
bins(ismember({kinData.mouse}, miceToExclude)) = 0; 
binLabels = conditions;
plotObsHeightTrajectories(kinData, bins.*includeTrial, binLabels, ['sensory cortex ' manipulation])
saveas(gcf, fullfile(getenv('OBSDATADIR'), 'figures', manipulation, 'SenKinematics.png'));
savefig(fullfile(getenv('OBSDATADIR'), 'figures', manipulation, 'SenKinematics.fig'))


