function sessions = getAllExperimentSessions()

% returns a cell array containing all sessions listed in the
% experimentMetaData spreadsheet. use this to get a list of all sessions if
% some variable needs to be recomputed for every session, for example

% settings
spreadsheet = fullfile(getenv('OBSDATADIR'), 'spreadSheets', 'experimentMetadata.xlsx');
sheets = {'baselineNotes', 'lesionNotes', 'senLesionNotes', 'muscimolNotes', 'sensoryDependenceNotes', 'whiskerTrimNotes'};  % which sheets to analyze

sessions = cell(1,length(sheets));
for i = 1:length(sheets)
    sheet = readtable(spreadsheet, 'Sheet', sheets{i});
    sessions{i} = sheet.session(~cellfun(@isempty, sheet.session));
end

sessions = unique(cat(1, sessions{:}));
fprintf('found %i sessions\n', length(sessions));