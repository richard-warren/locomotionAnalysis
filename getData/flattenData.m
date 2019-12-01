function dataOut = flattenData(data, varsToGet, varsGotten)

% todo: don't go deeper in tree if all vars are already found // how to
% deal with bins of continuous variables? // add conditionals here?

if ~exist('varsGotten', 'var'); varsGotten = struct(); end
if ischar(varsToGet); varsToGet={varsToGet}; end  % allows user to pass a single variable as a string, rather multiple variables using a cell array of strings


% determine which fields are structures
fields = fieldnames(data);
isFieldStruct = false(1,length(fields));
for i = 1:length(fields); isFieldStruct(i) = isstruct(data(1).(fields{i})); end

if any(ismember(fields, varsToGet)) || any(isFieldStruct)

    % store vars
    dataOut = repmat(varsGotten, 1, length(data)); % make duplicates of the rows to accomodate new data
    for field = fields(~isFieldStruct & ismember(fields, varsToGet)')'
        [dataOut(1:length(data)).(field{1})] = data.(field{1});
    end

    % get vars within nested structures
    if any(isFieldStruct)
        structData = cell(1,length(data));
        for row = 1:length(data)
            for field = fields(isFieldStruct)'
                structData{row} = flattenData(data(row).(field{1}), varsToGet, dataOut(row));
            end
        end
        dataOut = [structData{:}];
    end
    
else
    dataOut = varsGotten;
end

