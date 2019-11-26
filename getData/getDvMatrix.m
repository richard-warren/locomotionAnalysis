function dvMatrix = getDvMatrix(data, dv, vars, varsToAvg, conditionals, conditionInds)


% todo: multiple dvs // how to bin variables // return values of fields that are avged over, eg mouse name // meaningful error messages... // document, explain code logic

% analysisType is either 'mean' or 'median'


% initializations
if ~exist('conditionals', 'var'); conditionals = struct('name', ''); end
if ~exist('conditionInds', 'var'); conditionInds = nan(1, length(vars)); end
fields = fieldnames(data);
dvFound = ismember(dv, fields);

% initialize dvMatrix if dv is found
if ismember(dv, fields)
    dvMatrixSize = num2cell([cellfun(@length, {vars.levels}) length(data)]);
    dvMatrix = nan(dvMatrixSize{:});

% otherwise, find fields containing nested structs to be searched
else
    dvMatrices = cell(1,length(data));
    structFieldInds = false(1,length(fields));
    for i = 1:length(fields); structFieldInds(i) = isstruct(data(1).(fields{i})); end
    structFieldInds = find(structFieldInds);
end

% get inds of vars in data
varInds = find(ismember({vars.name}, fields));
avgDvs = any(ismember(varsToAvg, fields));
conditionalInds = find(ismember({conditionals.name}, fields));
conditionIndsSub = conditionInds;



% return empty matrix if branch contains neither dv nor nested structs
if ~dvFound && isempty(structFieldInds)
    dvMatrix = [];

% otherwise, loop through rows
else
    % find rows that meet conditionals
    binsToAnalyze = true(1,length(data));
    if ~isempty(conditionalInds)    
        for i = conditionalInds
            if any(cellfun(@ischar, ({data.(conditionals(i).name)}))) % if row contains strings
                binsToAnalyze = binsToAnalyze & ...
                    feval(conditionals(i).condition, {data.(conditionals(i).name)});
            else
                binsToAnalyze = binsToAnalyze & ...
                    feval(conditionals(i).condition, [data.(conditionals(i).name)]);
            end
        end
    end
    
    
    for i = find(binsToAnalyze)
        
%         try; disp(data(i).session); catch; end
%         try; if strcmp(data(i).mouse, 'sen11'); keyboard; end; catch; end

        % add vars to conditionInds
        skipRow = false;
        for j = varInds
            rowConditionInd = find(ismember(vars(j).levels, data(i).(vars(j).name)));
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
%                     try; if strcmp(data(i).session, '190326_001'); keyboard; end; catch; end
                    dvMatrices{i} = getDvMatrix(data(i).(fields{j}), ...
                        dv, vars, varsToAvg, conditionals, conditionIndsSub);
                    if ~isempty(dvMatrices{i}); break; end % terminate loop once a branch containing dv is found
                end

                if avgDvs
                    dvMatrices{i} = nanmean(dvMatrices{i}, length(vars)+1);
                end % avg matrix if data contains field in varsToAvg
            end
        end
    end
    
    if ~dvFound
%         dvMatrices = dvMatrices(~cellfun(@isempty, dvMatrices));  %  remove empty entries (not sure if it is safe to do this)
        dvMatrix = cat(length(vars)+1, dvMatrices{:});
    end % concatenate matrices obtained from nested structs
end



