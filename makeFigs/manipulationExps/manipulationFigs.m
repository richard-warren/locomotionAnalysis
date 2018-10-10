%% load session info


% settings
manipulation = 'lesion';
maxLesionSession = 4;


if strcmp(manipulation, 'muscimol'); conditions = {'saline', 'muscimol'};
elseif strcmp(manipulation, 'lesion'); conditions = {'pre', 'post'}; end
sessionInfo = readtable([getenv('OBSDATADIR') 'sessions\sessionInfo.xlsx'], 'Sheet', [manipulation 'Notes']);
sessionInfo = sessionInfo(~cellfun(@isempty, sessionInfo.session) & strcmp(sessionInfo.brainRegion, 'mtc'),:);
sessionInfo = sessionInfo(~cellfun(@isempty, sessionInfo.session) & logical([sessionInfo.include]),:);



%% compute kinematic data

loadPreviousData = false;

fileName = [getenv('OBSDATADIR') 'matlabData\' manipulation 'kinematicData.mat'];
if loadPreviousData && exist(fileName, 'file')
    load(fileName, 'data');
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


%%





% explore relationship between speed, success, and manipulation






%% plot histograms
controlBins = strcmp({speedAvoidanceData.condition}, conditions{1});
manipBins = strcmp({speedAvoidanceData.condition}, conditions{2});

figure;
histogram([speedAvoidanceData(controlBins).avgVel], 50); hold on
histogram([speedAvoidanceData(manipBins).avgVel], 50);
legend(conditions)

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

%% find matched bins

tolerance = .02;

binEdges = 0:tolerance:max([speedAvoidanceData.avgVel]);
[counts, binEdges, bins] = histcounts([speedAvoidanceData.avgVel], binEdges);


if strcmp(manipulation, 'lesion'); includeTrial = [speedAvoidanceData.conditionNum]<=maxLesionSession;
else; includeTrial = true(1,length(speedAvoidanceData)); end

controlBins = strcmp({speedAvoidanceData.condition}, conditions{1});
manipBins = strcmp({speedAvoidanceData.condition}, conditions{2});
matchedBins = false(1,length(speedAvoidanceData));
mice = unique({speedAvoidanceData.mouse});

for i = 1:length(binEdges-1)
    for j = 1:length(mice)
    
        mouseBins = strcmp({speedAvoidanceData.mouse}, mice{j});
        binControlTrials = find(bins==i & controlBins & mouseBins & includeTrial);
        binManipTrials = find(bins==i & manipBins & mouseBins & includeTrial);

        if ~isempty(binControlTrials) && ~isempty(binManipTrials)
            minCount = min(length(binControlTrials), length(binManipTrials));
            matchedBins(binControlTrials(randperm(length(binControlTrials), minCount))) = true;
            matchedBins(binManipTrials(randperm(length(binManipTrials), minCount))) = true;
        end
    end
end


% plot histograms
figure;
histogram([speedAvoidanceData(controlBins & matchedBins).avgVel], 50); hold on
histogram([speedAvoidanceData(manipBins & matchedBins).avgVel], 50);
legend(conditions)












