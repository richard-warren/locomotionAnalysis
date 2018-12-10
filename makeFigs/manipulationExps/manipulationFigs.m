%% load session info


% settings
manipulation = 'lesion';
maxLesionSession = 3;


if strcmp(manipulation, 'muscimol'); conditions = {'saline', 'muscimol'};
elseif strcmp(manipulation, 'lesion'); conditions = {'pre', 'post', 'postNoWisk'}; end
sessionInfo = readtable([getenv('OBSDATADIR') 'sessions\sessionInfo.xlsx'], 'Sheet', [manipulation 'Notes']);
% sessionInfo = sessionInfo(~cellfun(@isempty, sessionInfo.session) & strcmp(sessionInfo.brainRegion, 'mtc'),:);
% sessionInfo = sessionInfo(~cellfun(@isempty, sessionInfo.session) & strcmp(sessionInfo.brainRegion, 'sen'),:);
sessionInfo = sessionInfo(~cellfun(@isempty, sessionInfo.session) & [sessionInfo.include]==1,:); % remove empty entires in spreadsheet



%% compute kinematic data
tic
loadPreviousData = false;

fileName = [getenv('OBSDATADIR') 'matlabData\' manipulation 'kinematicData.mat'];
if loadPreviousData && exist(fileName, 'file')
    load(fileName, 'data');
    kinData = getKinematicData4(sessionInfo.session, sessionInfo, data);
else
    kinData = getKinematicData4(sessionInfo.session, sessionInfo, []);
end
data = kinData; save([getenv('OBSDATADIR') 'matlabData\' manipulation 'KinematicData.mat'], 'data', '-v7.3', '-nocompression'); clear data;

if strcmp(manipulation, 'lesion'); kinData = kinData([kinData.conditionNum]<=maxLesionSession); end
toc

%% load kinematic data

load([getenv('OBSDATADIR') 'matlabData\' manipulation 'KinematicData.mat'], 'data');
kinData = data; clear data;
if strcmp(manipulation, 'lesion'); kinData = kinData([kinData.conditionNum]<=maxLesionSession); end
disp([manipulation ' kinematic data loaded!'])

%% compute speed and avoidance data

speedAvoidanceData = getSpeedAndObsAvoidanceData(sessionInfo.session, sessionInfo, false);
data = speedAvoidanceData; save([getenv('OBSDATADIR') 'matlabData\' manipulation 'SpeedAvoidanceData.mat'], 'data'); clear data;
if strcmp(manipulation, 'lesion'); kinData = speedAvoidanceData([speedAvoidanceData.conditionNum]<=maxLesionSession); end


%% load speed and avoidance data

load([getenv('OBSDATADIR') 'matlabData\' manipulation 'SpeedAvoidanceData.mat'], 'data');
speedAvoidanceData = data; clear data;
if strcmp(manipulation, 'lesion'); kinData = speedAvoidanceData([speedAvoidanceData.conditionNum]<=maxLesionSession); end
disp([manipulation ' speed avoidance data loaded!'])

%% find matched bins

% settings
velTolerance = .05;
angleTolerance = 2;

[matchedBinsSpeedAvoidanceBins, weightsSpeedAvoidance] = findMatchedBins(speedAvoidanceData, conditions, velTolerance, angleTolerance);


% % plot matched bins histogram
% figure;
% angle = abs([speedAvoidanceData.avgAngle]);
% speed = [speedAvoidanceData.avgVel];
% hist3([angle(controlBins & matchedBins)', speed(controlBins & matchedBins)'], 'FaceColor', [1 0 0], 'FaceAlpha', 0.5); hold on
% hist3([angle(manipBins & matchedBins)', speed(manipBins & matchedBins)'], 'FaceColor', [0 1 0], 'FaceAlpha', 0.5);



%% ----------
% PLOT THINGS
%  ----------

%% plot dv averages

matchDistributions = true;


if matchDistributions; weights = weightsSpeedAvoidance; else; weights = ones(1,length(speedAvoidanceData)); end

manipulationBarPlots(speedAvoidanceData, conditions, manipulation, weights);
saveas(gcf, [getenv('OBSDATADIR') 'figures\' manipulation '\' manipulation 'BarPlots.png']);
savefig([getenv('OBSDATADIR') 'figures\' manipulation '\' manipulation 'BarPlots.fig'])


%% plot dvs across sessions

manipulationAcrossSessions(speedAvoidanceData, conditions, manipulation);
saveas(gcf, [getenv('OBSDATADIR') 'figures\' manipulation '\' manipulation 'AcrossSession.png']);
savefig([getenv('OBSDATADIR') 'figures\' manipulation '\' manipulation 'AcrossSessions.fig'])

%% paw height by obs height for mtc

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


%% paw height by obs height for sen

binNames = conditions;

velBins = zeros(1,length(kinData));
velBins(strcmp({kinData.condition}, conditions{1}) & strcmp({kinData.brainRegion}, 'sen')) = 1; % pre
velBins(strcmp({kinData.condition}, conditions{2}) & strcmp({kinData.brainRegion}, 'sen')) = 2; % post
if length(conditions)==3; velBins(strcmp({kinData.condition}, conditions{3}) & strcmp({kinData.brainRegion}, 'sen')) = 3; end % post no wisk

% don't include too many post lesion sessions if manipulation=='lesion'
if strcmp(manipulation, 'lesion'); includeTrial = [kinData.conditionNum]<=maxLesionSession;
else; includeTrial = ones(1,length(kinData)); end

scatterObsVsPawHeights(kinData, velBins.*includeTrial, binNames);
saveas(gcf, fullfile(getenv('OBSDATADIR'), 'figures/', manipulation, '/heightShapingScatterSen.png'));
savefig(fullfile(getenv('OBSDATADIR'), 'figures/', manipulation, '/heightShapingScatterSend.fig'))


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

%% sen
velBins = zeros(1,length(kinData));
velBins(senBins & controlBins) = 1; % control
velBins(senBins & manipBins) = 2; % manipluation
velBins(senBins & noWiskBins) = 3; % post no wisk
velBins(ismember({kinData.mouse}, miceToExclude)) = 0; 
binLabels = conditions;
plotObsHeightTrajectories(kinData, velBins.*includeTrial, binLabels, ['sensory cortex ' manipulation])
saveas(gcf, fullfile(getenv('OBSDATADIR'), 'figures', manipulation, 'SenKinematics.png'));
savefig(fullfile(getenv('OBSDATADIR'), 'figures', manipulation, 'SenKinematics.fig'))


%%





% explore relationship between speed, success, and manipulation






%% plot histograms
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





