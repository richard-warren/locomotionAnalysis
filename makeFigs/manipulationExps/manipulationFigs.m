%% load session info


% settings
manipulation = 'muscimol';
brainRegion = 'mtc';
maxLesionSession = 3;


% load session metadata
sessionInfo = readtable(fullfile(getenv('OBSDATADIR'), 'spreadSheets', 'experimentMetadata.xlsx'), 'Sheet', [manipulation 'Notes']);
sessionInfo = sessionInfo(~cellfun(@isempty, sessionInfo.session) & ...
                          [sessionInfo.include]==1,:); % remove empty rows and sessions not to be included
sessionInfo = sessionInfo(strcmp(sessionInfo.brainRegion, brainRegion),:); % keep data for single brain region only

%% compute kinData for all sessions
sessions = unique(sessionInfo.session);
parfor i = 1:length(sessions)
    getKinematicData5(sessions{i});
end

%% compute kinematic data
loadPreviousData = false;

fileName = fullfile(getenv('OBSDATADIR'), 'matlabData', [brainRegion '_' manipulation '_kinematicData.mat']);
if loadPreviousData && exist(fileName, 'file')
    load(fileName, 'data');
    kinData = getKinematicData4(sessionInfo.session, sessionInfo, data);
else
    kinData = getKinematicData4(sessionInfo.session, sessionInfo, []);
end
data = kinData; save(fileName, 'data', '-v7.3', '-nocompression'); clear data;
% if strcmp(manipulation, 'lesion'); kinBins=[kinData.conditionNum]<=maxLesionSession; else; bins=true(size(kinData)); end

%% load kinematic data

load(fullfile(getenv('OBSDATADIR'), 'matlabData', [brainRegion '_' manipulation '_kinematicData.mat']), 'data');
kinData = data; clear data;
% if strcmp(manipulation, 'lesion'); kinBins=[kinData.conditionNum]<=maxLesionSession; else; bins=true(size(kinData)); end
disp([manipulation ' kinematic data loaded!'])

%% compute speed and avoidance data

speedAvoidanceData = getSpeedAndObsAvoidanceData(sessionInfo.session, sessionInfo, true);
data = speedAvoidanceData; save(fullfile(getenv('OBSDATADIR'), 'matlabData', [brainRegion '_' manipulation '_speedAvoidanceData.mat']), 'data'); clear data;
% if strcmp(manipulation, 'lesion'); speedBins=[speedAvoidanceData.conditionNum]<=maxLesionSession; else; bins=true(size(speedAvoidanceData)); end

%% load speed and avoidance data

load(fullfile(getenv('OBSDATADIR'), 'matlabData', [brainRegion '_' manipulation '_speedAvoidanceData.mat']), 'data');
speedAvoidanceData = data; clear data;
% if strcmp(manipulation, 'lesion'); speedBins=[speedAvoidanceData.conditionNum]<=maxLesionSession; else; bins=true(size(speedAvoidanceData)); end
disp([manipulation ' speed avoidance data loaded!'])

%% find matched bins

% settings
velTolerance = .05;
angleTolerance = 2;


[matchedBinsSpeedAvoidanceBins, weightsSpeedAvoidance] = findMatchedBins(speedAvoidanceData(speedBins), conditions, velTolerance, angleTolerance);


% % plot matched bins histogram
% figure;
% angle = abs([speedAvoidanceData.avgAngle]);
% speed = [speedAvoidanceData.avgVel];
% hist3([angle(controlBins & matchedBins)', speed(controlBins & matchedBins)'], 'FaceColor', [1 0 0], 'FaceAlpha', 0.5); hold on
% hist3([angle(manipBins & matchedBins)', speed(manipBins & matchedBins)'], 'FaceColor', [0 1 0], 'FaceAlpha', 0.5);

%%  compute dependent measures

dvs = {'success', 'forePawErrRate', 'speed', 'bodyAngleContra', 'contraFirstRate', 'bigStepProb', 'pawHgt', 'hgtShaping'};

sessionDvs = getSessionDvs(dvs, speedAvoidanceData, kinData);
disp('finished computing dependent measures!')

%% ----------
% PLOT THINGS
%  ----------

%% bar plots

dvs = {'success', 'speed', 'bodyAngleContra', 'forePawErrRateIpsi', 'forePawErrRateContra', ...
    'contraFirstRate', 'bigStepProbIpsi', 'bigStepProbContra', 'pawHgtIpsi', 'pawHgtContra', ...
    'hgtShapingIpsi', 'hgtShapingContra'};

if strcmp(manipulation, 'lesion'); bins = [sessionDvs.conditionNum]<=maxLesionSession; else; bins = true(size(sessionDvs)); end
barPlots(sessionDvs(bins), dvs, [brainRegion '_' manipulation], flipud(conditions))
saveas(gcf, fullfile(getenv('OBSDATADIR'), 'figures', 'manipulations', [brainRegion '_' manipulation 'BarPlots.png']));
savefig(fullfile(getenv('OBSDATADIR'), 'figures', 'manipulations', [brainRegion '_' manipulation 'BarPlots.fig']));


%% sessions over time plots

miceToShow = 'all'; % set to 'all' to show all mice

if strcmp(miceToShow, 'all'); bins = true(1,length(sessionDvs)); else; bins = ismember({sessionDvs.mouse}, miceToShow); end
plotAcrossSessions(sessionDvs, dvs, [brainRegion '_' manipulation])
saveas(gcf, fullfile(getenv('OBSDATADIR'), 'figures', 'manipulations', [brainRegion '_' manipulation 'AcrossSessions.png']));
savefig(fullfile(getenv('OBSDATADIR'), 'figures', 'manipulations', [brainRegion '_' manipulation 'AcrossSessions.fig']));

%% speed vs. position plots

if strcmp(manipulation, 'lesion'); bins = [speedAvoidanceData.conditionNum]<=maxLesionSession; else; bins = true(size(speedAvoidanceData)); end
plotSpeedVsPosition(speedAvoidanceData(bins), [brainRegion '_' manipulation])
saveas(gcf, fullfile(getenv('OBSDATADIR'), 'figures', 'manipulations', [brainRegion '_' manipulation 'SpeedVsPos.png']));
savefig(fullfile(getenv('OBSDATADIR'), 'figures', 'manipulations', [brainRegion '_' manipulation 'SpeedVsPos.fig']));

%% speed at wisk contact

plotSpeedAtWiskContact(sessionInfo(sessionInfo.conditionNum<=maxLesionSession,:), [brainRegion '_' manipulation])
saveas(gcf, fullfile(getenv('OBSDATADIR'), 'figures', 'manipulations', [brainRegion '_' manipulation 'SpeedAtWiskContact.png']));
savefig(fullfile(getenv('OBSDATADIR'), 'figures', 'manipulations', [brainRegion '_' manipulation 'SpeedAtWiskContact.fig']));


%% !!! paw height by obs height

binNames = {'ipsi', 'contra', conditions{1}};

velBins = zeros(1,length(kinData));
velBins(strcmp({kinData.condition}, conditions{2}) & strcmp({kinData.brainRegion}, 'mtc') & [kinData.ipsiPawFirst]) = 1;
velBins(strcmp({kinData.condition}, conditions{2}) & strcmp({kinData.brainRegion}, 'mtc') & [kinData.contraPawFirst]) = 2;
velBins(strcmp({kinData.condition}, conditions{1}) & strcmp({kinData.brainRegion}, 'mtc')) = 3;

% don't include too many post lesion sessions if manipulation=='lesion'
if strcmp(manipulation, 'lesion'); includeTrial = [kinData.conditionNum]<=maxLesionSession;
else; includeTrial = ones(1,length(kinData)); end

scatterObsVsPawHeights(kinData, velBins.*includeTrial, binNames);
saveas(gcf, fullfile(getenv('OBSDATADIR'), 'figures/', manipulation, '/heightShapingScatterMtc.png'));
savefig(fullfile(getenv('OBSDATADIR'), 'figures/', manipulation, '/heightShapingScatterMtc.fig'))


%% plot kinematics

% settings
% miceToExclude = {'sen5'};
miceToExclude = {''};

% get trial bins
contraFirstBins = [kinData.contraPawFirst];
ipsiFirstBins = [kinData.ipsiPawFirst];
controlBins = strcmp({kinData.condition}, conditions{1});
manipBins = strcmp({kinData.condition}, conditions{2});
noWiskBins = strcmp({kinData.condition}, 'postNoWisk');
lightOnBins = [kinData.isLightOn];
mtcBins = strcmp({kinData.brainRegion}, 'mtc');
senBins = strcmp({kinData.brainRegion}, 'sen');

% don't include too many post lesion sessions if manipulation=='lesion'
if strcmp(manipulation, 'lesion'); includeTrial = [kinData.conditionNum]<=maxLesionSession;
else; includeTrial = ones(1,length(kinData)); end

% mtc
velBins = zeros(1,length(kinData));
velBins(mtcBins & ipsiFirstBins & manipBins) = 1; % ipsi
velBins(mtcBins & contraFirstBins & manipBins) = 2; % contra
velBins(mtcBins & controlBins) = 3; % pre
binLabels = {'ipsi', 'contra', 'pre'};
plotObsHeightTrajectories(kinData, velBins.*includeTrial, binLabels, ['motor cortex ' manipulation])
saveas(gcf, fullfile(getenv('OBSDATADIR'), 'figures', manipulation, 'MtcKinematics.png'));
savefig(fullfile(getenv('OBSDATADIR'), 'figures', manipulation, 'MtcKinematics.fig'))







%% explore relationship between speed, success, and manipulation

% plot histograms
controlBins = strcmp({speedAvoidanceData.condition}, conditions{1});
manipBins = strcmp({speedAvoidanceData.condition}, conditions{2});

figure;
histogram([speedAvoidanceData(controlBins).avgVel], 50); hold on
histogram([speedAvoidanceData(manipBins).avgVel], 50);
legend(conditions)

figure;
histogram(abs([speedAvoidanceData(controlBins).avgAngle]), 50); hold on
histogram(abs([speedAvoidanceData(manipBins).avgAngle]), 50);
legend(conditions)

figure;
angle = abs([speedAvoidanceData.avgAngle]);
speed = [speedAvoidanceData.avgVel];
hist3([angle(controlBins)', speed(controlBins)'], 'FaceColor', [1 0 0], 'FaceAlpha', 0.5); hold on
hist3([angle(manipBins)', speed(manipBins)'], 'FaceColor', [0 1 0], 'FaceAlpha', 0.5);




%% plot succcses vs speed

dv = .01;
width = .05;
velLims = [0 1];
successThresh = 5;
validBins = strcmp({speedAvoidanceData.condition}, conditions{1}); % only look at control condition trials

isSuccess = cellfun(@sum, {speedAvoidanceData.totalTouchFramesPerPaw}) < successThresh;
binCenters = velLims(1)+.5*width : dv : velLims(2)-.5*width;

successRates = nan(1, length(binCenters));

for i = 1:length(binCenters)
    
    trialBins = (validBins & ...
                 [speedAvoidanceData.avgVel] > binCenters(i)-.5*width & ...
                 [speedAvoidanceData.avgVel] < binCenters(i)+.5*width);
    if any(trialBins)
        successRates(i) = mean(isSuccess(trialBins));
    end
end

figure; plot(binCenters, successRates)





