%% load session info

sessionInfo = readtable([getenv('OBSDATADIR') 'sessions\sessionInfo.xlsx'], 'Sheet', 'baselineNotes');
sessionInfo = sessionInfo(sessionInfo.include==1 & ~cellfun(@isempty, sessionInfo.session),:);



%% get kinematic data
obsPos = -0.0087;
loadPreviousData = false;

if loadPreviousData
    load([getenv('OBSDATADIR') 'matlabData\baselineKinematicData.mat'], 'data');
    kinData = getKinematicData4(sessionInfo.session, data, obsPos);
else
    kinData = getKinematicData4(sessionInfo.session, [], obsPos);
end
data = kinData; save([getenv('OBSDATADIR') 'matlabData\baselineKinematicData.mat'], 'data');
