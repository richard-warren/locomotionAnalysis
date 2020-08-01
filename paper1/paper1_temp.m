% re-analyze no whisker sessions

noWiskExperiments = {'senLesionNotes', 'sensoryDependenceNotes'};

for i = 1:length(noWiskExperiments)
    
    expData = readtable(fullfile(getenv('OBSDATADIR'), 'spreadSheets', 'experimentMetaData.xlsx'), 'Sheet', noWiskExperiments{i});
    expData = expData(expData.include==1,:);
    noWiskSessions = expData.session(strcmp(expData.whiskers, 'none'));
    
    for j = 1:length(noWiskSessions)
        fprintf('\n-------- %s --------\n', noWiskSessions{j})
        analyzeSession(noWiskSessions{j}, 'overwriteVars', 'wiskContactFrames');
        getKinematicData(noWiskSessions{j});
    end
end