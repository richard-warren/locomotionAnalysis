function makeCellDataCsv(session, varargin)

s.fileName = 'cellData.csv';


% initializations
if exist('varargin', 'var'); for i = 1:2:length(varargin); s.(varargin{i}) = varargin{i+1}; end; end  % parse name-value pairs
fullFileName = fullfile(getenv('OBSDATADIR'), 'sessions', session, s.fileName);

% check if these unit_ids have already been written to cellData.csv file
if exist(fullFileName, 'file')
    fprintf('WARNING! %s already exists. %s not created.\n', ...
        fullFileName, s.fileName)
else
    fprintf('creating %s\n', fullFileName)
    csvTable = cell2table(cell(0,6), ...
        'VariableNames', {'unit_id', 'include', 'timeStart', 'timeEnd', 'location', 'notes'});
    [~, unit_ids] = getGoodSpkInds(session);
    warning('off', 'MATLAB:table:RowsAddedExistingVars')
    csvTable(1:length(unit_ids),1) = table(unit_ids);
    writetable(csvTable, fullFileName)
    warning('on', 'MATLAB:table:RowsAddedExistingVars')
end
