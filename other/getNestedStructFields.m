function dataOut = getNestedStructFields(data, varsToGet, varsGotten)

% to do: don't go deeper in tree if all vars are already gotten

if ~exist('varsGotten', 'var'); varsGotten = struct(); end


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
                structData{row} = getNestedStructFields(data(row).(field{1}), varsToGet, dataOut(row)); %catch; keyboard; end
            end
        end
        dataOut = [structData{:}];
    end
    
else
    dataOut = varsGotten;
end

