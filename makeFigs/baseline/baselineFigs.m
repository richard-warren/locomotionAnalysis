%% load session info

sessionInfo = readtable([getenv('OBSDATADIR') 'sessions\sessionInfo.xlsx'], 'Sheet', 'baselineNotes');
sessionInfo = sessionInfo(sessionInfo.include==1 & ~cellfun(@isempty, sessionInfo.session),:);



%% compute kinematic data

loadPreviousData = false;

if loadPreviousData
    load([getenv('OBSDATADIR') 'matlabData\baselineKinematicData.mat'], 'data');
    kinData = getKinematicData4(sessionInfo.session, sessionInfo, data);
else
    kinData = getKinematicData4(sessionInfo.session, sessionInfo, []);
end
data = kinData; save([getenv('OBSDATADIR') 'matlabData\baselineKinematicData.mat'], 'data');

%% load previous data

load([getenv('OBSDATADIR') 'matlabData\baselineKinematicData.mat'], 'data');
kinData = data; clear data;

%% plot one step prob by obs height

oneStepProbByObsHeight(kinData);
saveas(gcf, fullfile(getenv('OBSDATADIR'), 'figures/baseline/oneStepProbByObsHeight.png'));
savefig(fullfile(getenv('OBSDATADIR'), 'figures/baseline/oneStepProbByObsHeight.fig'))

%% plot one vs two step trajectories

binNum = 3;
binVar = [kinData.swingStartDistance] + [kinData.predictedLengths]; % predicted distance to obs
binVar(abs(zscore(binVar))>3) = 0; % remove outliers
binEdges = linspace(min(binVar), max(binVar), binNum+1);
bins = discretize(binVar, binEdges);
binLabels = cell(1,binNum);
for i = 1:binNum; binLabels{i} = sprintf('%.3f', mean(binVar(bins==i))); end

plotOneVsTwoStepTrajectories(kinData, bins, 'averages')
saveas(gcf, fullfile(getenv('OBSDATADIR'), 'figures/baseline/oneVsTwoStepTrajectories.png'));
savefig(fullfile(getenv('OBSDATADIR'), 'figures/baseline/oneVsTwoStepTrajectories.fig'))

%% plot obs height kinematics

binNum = 1;
binVar = [kinData.swingStartDistance] + [kinData.predictedLengths]; % predicted distance to obs
binVar(abs(zscore(binVar))>3) = 0; % remove outliers
binEdges = linspace(min(binVar), max(binVar), binNum+1);
bins = discretize(binVar, binEdges);
binLabels = cell(1,binNum);
for i = 1:binNum; binLabels{i} = sprintf('%.3f', mean(binVar(bins==i))); end

if binNum==1; binLabels={''}; end
plotObsHeightTrajectories(kinData, bins, binLabels, 'baselineHeightKinematics')
saveas(gcf, fullfile(getenv('OBSDATADIR'), 'figures/baseline/heightKinematics.png'));
savefig(fullfile(Wgetenv('OBSDATADIR'), 'figures/baseline/heightKinematics.fig'))

%% scatter obs vs paw heights

scatterObsVsPawHeights(kinData, ones(1,length(kinData)), {''})
saveas(gcf, fullfile(getenv('OBSDATADIR'), 'figures/baseline/obsVsPawHeights.png'));
savefig(fullfile(getenv('OBSDATADIR'), 'figures/baseline/obsVsPawHeights.fig'))






