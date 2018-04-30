function sessions = selectSessions

sessionDirs = uigetdir2([getenv('OBSDATADIR') 'sessions\'], 'select folders to analyze');
sessions = cellfun(@(x) x(end-9:end), sessionDirs, 'uniformoutput', false);
sessions = strjoin(sessions, ' ');
disp(sessions)