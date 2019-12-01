function [sessions, experiments] = getAllExperimentSessions(sheets)

% returns a cell array containing all sessions listed in the
% experimentMetaData spreadsheet. use this to get a list of all sessions if
% some variable needs to be recomputed for every session, for example. also
% retruns experiments, which lists the names of the sheet from which the
% session was drawn

% settings
spreadsheet = fullfile(getenv('OBSDATADIR'), 'spreadSheets', 'experimentMetadata.xlsx');
if ~exist('sheets', 'var')
    sheets = {'baselineNotes', 'mtcLesionNotes', 'senLesionNotes', 'muscimolNotes', 'sensoryDependenceNotes', 'whiskerTrimNotes'};  % which sheets to analyze
elseif ischar(sheets)
    sheets = {sheets};
end

sessions = cell(1,length(sheets));
experiments = cell(1,length(sheets));

for i = 1:length(sheets)
    sheet = readtable(spreadsheet, 'Sheet', sheets{i});
    bins = ~cellfun(@isempty, sheet.session);
    sessions{i} = sheet.session(bins);
    experiments{i} = repelem(sheets(i), sum(bins))';
end

% keyboard
[sessions, inds] = unique(cat(1, sessions{:}));
experiments = cat(1, experiments{:});
experiments = {experiments{inds}};
fprintf('found %i sessions\n', length(sessions));