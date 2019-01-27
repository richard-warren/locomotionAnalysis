function dvMatrix = barPlotWrapper(data, dv, vars, varLevels, varLevelNames, varsToAvg, conditionInds)


% TO DO: nest vars, varLevels, and varLevelNames? // meaningful error
% messages... // document, explain code logic


% initializations
if ~exist('conditionInds', 'var'); conditionInds = nan(1, length(vars)); end
fields = fieldnames(data);
dvFound = ismember(dv, fields);

% initialize dvMatrix if dv is found
if ismember(dv, fields)
    dvMatrixSize = num2cell([cellfun(@length, varLevels) length(data)]);
    dvMatrix = nan(dvMatrixSize{:});

% otherwise, find fields containing nested structs to be searched
else
    dvMatrices = cell(1,length(data));
    structFieldInds = false(1,length(fields));
    for i = 1:length(fields); structFieldInds(i) = isstruct(data(1).(fields{i})); end
    structFieldInds = find(structFieldInds);
end

% get inds of vars in data
varInds = find(ismember(vars, fields));
avgDvs = any(ismember(varsToAvg, fields));
conditionIndsSub = conditionInds;



% return empty matrix if branch contains neither dv nor nested structs
if ~dvFound && isempty(structFieldInds)
    dvMatrix = [];

% otherwise, loop through rows
else
    for i = 1:length(data)

        % add vars to conditionInds
        skipRow = false;
        for j = varInds
            rowConditionInd = find(ismember(varLevels{j}, data(i).(vars{j})));
            if ~isempty(rowConditionInd)
                conditionIndsSub(j) = rowConditionInd;
            else
                skipRow = true;
            end
        end
        
        
        if ~skipRow
            % if dv is in data, add to correct location in dvMatrix
            if dvFound
                dvMatrixInds = num2cell([conditionIndsSub i]);
                dvMatrix(dvMatrixInds{:}) = data(i).(dv);

            % otherwise loop through nested structs
            else
                for j = structFieldInds
                    dvMatrices{i} = barPlotWrapper(data(i).(fields{j}), ...
                        dv, vars, varLevels, varLevelNames, varsToAvg, conditionIndsSub);
                    if ~isempty(dvMatrices{i}); break; end % terminate loop once a branch containing dv is found
                end

                if avgDvs; dvMatrices{i} = nanmean(dvMatrices{i}, length(vars)+1); end % avg matrix if data contains field in varsToAvg
            end
        end
    end

    if ~dvFound; dvMatrix = cat(length(vars)+1, dvMatrices{:}); end % concatenate matrices obtained from nested structs
end



