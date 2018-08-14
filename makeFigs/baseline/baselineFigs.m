%% load session info

sessionInfo = readtable([getenv('OBSDATADIR') 'sessions\sessionInfo.xlsx'], 'Sheet', 'baselineNotes');
sessionInfo = sessionInfo(sessionInfo.include==1 & ~cellfun(@isempty, sessionInfo.session),:);



%% compute kinematic data
obsPos = -0.0087;
loadPreviousData = true;

if loadPreviousData
    load([getenv('OBSDATADIR') 'matlabData\baselineKinematicData.mat'], 'data');
    kinData = getKinematicData4(sessionInfo.session, data, obsPos);
else
    kinData = getKinematicData4(sessionInfo.session, [], obsPos);
end
data = kinData; save([getenv('OBSDATADIR') 'matlabData\baselineKinematicData.mat'], 'data');

%% load previous data

load([getenv('OBSDATADIR') 'matlabData\baselineKinematicData.mat'], 'data');
kinData = data; clear data;

%% plot one step prob by obs height

oneStepProbByObsHeight(kinData);
saveas(gcf, fullfile(getenv('OBSDATADIR'), 'figures/baseline/oneStepProbByObsHeight.png'));
savefig(fullfile(getenv('OBSDATADIR'), 'figures/baseline/oneStepProbByObsHeight.fig'))
