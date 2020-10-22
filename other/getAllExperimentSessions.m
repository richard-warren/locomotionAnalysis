function experimentSessions = getAllExperimentSessions(varargin)

% returns table of all sessions and mice in all experiments listed in
% s.experiments // reads this information from experimentMetadata.xlsx

% settings
s.experiments = {'baselineNotes', 'mtcLesionNotes', 'senLesionNotes', 'muscimolNotes', 'sensoryDependenceNotes', 'whiskerTrimNotes'};
s.includeOnly = true;  % whether to only include sessions where include column is true in spreadsheet
s.extraColumns = {};   % copy these columns for experiment spreadsheet into experimentSessions output



% initializations
if exist('varargin', 'var'); for i = 1:2:length(varargin); s.(varargin{i}) = varargin{i+1}; end; end  % reassign settings passed in varargin
columns = [{'mouse', 'session', 'include'} s.extraColumns];  % copy these fields into a composite table across experiments
spreadsheet = fullfile(getenv('OBSDATADIR'), 'spreadSheets', 'experimentMetadata.xlsx');
if ischar(s.experiments); s.experiments = {s.experiments}; end

experimentSessions = cell(1, length(s.experiments));
for i = 1:length(s.experiments)
    sheet = readtable(spreadsheet, 'Sheet', s.experiments{i});
    rowBins = ~cellfun(@isempty, sheet.session);
    [~, colInds] = ismember(columns, sheet.Properties.VariableNames);
    experimentColumn = table(repelem({s.experiments{i}}, sum(rowBins), 1), 'VariableNames', {'experiment'});
    experimentSessions{i} = [experimentColumn sheet(rowBins, colInds)];
end

experimentSessions = cat(1, experimentSessions{:});
if s.includeOnly; experimentSessions = experimentSessions(experimentSessions.include==1,:); end

fprintf('found %i sessions\n', height(experimentSessions));