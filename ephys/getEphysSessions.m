function sessions = getEphysSessions()

% returns all sessions in ephysInfo.xlsx for which include==1

ephysInfo = readtable(fullfile(getenv('OBSDATADIR'), 'spreadSheets', 'ephysInfo.xlsx'));
sessions = ephysInfo.session(ephysInfo.include==1);