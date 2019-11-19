function [sessions, sessionText] = selectSessions

sessionDirs = uigetdir2([getenv('OBSDATADIR') 'sessions\'], 'select folders to analyze');
sessions = cellfun(@(x) x(end-9:end), sessionDirs, 'uniformoutput', false);
sessionsText = strjoin(sessions, ' ');
% disp(sessionsText)