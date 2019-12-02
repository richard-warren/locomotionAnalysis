function dataOut = flattenData(data, varsToGet, varsGotten)

% given nested 'data' struct (see 'getExperimentData'), returns flattened,
% non-nested, data structure with one row per entry at the level of the
% deepest variable in varsToGet. varsToGet is cell array of variables that
% should be included as columns in dataOut. e.g., if you varsToGet =
% {'mouse', 'trialVel'}, each row will correspond to a trial, whereas if
% varsToGet = {'mouse', 'isFore'}, each row will correspond to a PAW,
% becuase isFore is a paw level varialbe. don't use varsGotten, which is
% only used in the recursive calls of this function

% todo: don't go deeper in tree if all vars are already found // how to
% deal with bins of continuous variables? // add conditionals here?


if ~exist('varsGotten', 'var'); varsGotten = struct(); end
if ischar(varsToGet); varsToGet={varsToGet}; end  % allows user to pass a single variable as a string, rather multiple variables using a cell array of strings


% determine which fields are structures
fields = fieldnames(data);
isStruct = false(1,length(fields));
for i = 1:length(fields); isStruct(i) = isstruct(data(1).(fields{i})); end

if any(ismember(fields, varsToGet)) || any(isStruct)

    % store vars
    dataOut = repmat(varsGotten, 1, length(data));  % make duplicates of the rows to accomodate new data
    for field = fields(~isStruct & ismember(fields, varsToGet)')'
        [dataOut(1:length(data)).(field{1})] = data.(field{1});
    end

    % get vars within nested structures
    if any(isStruct)  && ~all(ismember(varsToGet, fieldnames(dataOut)))  % if there are any further nested structures, and we haven't gotten all of the variables already
        structData = cell(1,length(data));
        for row = 1:length(data)
            for field = fields(isStruct)'
                structData{row} = flattenData(data(row).(field{1}), varsToGet, dataOut(row));
            end
        end
        try dataOut = [structData{:}]; catch; keyboard; end  % todo: this will break if a variable is found in some nested structs, but not others... // i should find intersection between fielnames of all structs, and keep only columns in the intersection, and give error message saying that some varaibles were ommitted
    end
    
else
    dataOut = varsGotten;
end

