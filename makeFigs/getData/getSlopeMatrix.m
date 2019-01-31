function dvMatrix = getSlopeMatrix(data, dvs, vars, varsToAvg, varsToSlope, ...
                                   conditionals, conditionInds, dvsFound, isDvFound)

% given nested struct data, finds TWO dvs, and computes the slope between
% these dvs at the level of varsToSlope // eg, if varsToSlope is 'session',
% it gets all of the two DVs for a given session, then computes the slope
% for each session // dvMatrix has dimension [var1, var2, 2 (dv num), sample
% num] // after slopes are computed, second to last dimension has length 1,
% as the two vectors are replaced by a single number representing the slope
% between the two vectors... god this is a mess

% initializations
isFirstPass = ~exist('conditionInds', 'var');
if isFirstPass
    conditionInds = nan(1, length(vars));
    dvsFound = [nan nan];
    isDvFound = [false false];
    if ~exist('conditionals', 'var'); conditionals = struct('name', ''); end
end
fields = fieldnames(data);
dimNum = length(vars)+2; % one extra dim for the two dvs, and one for the individual samples


% find fields containing nested structs to be searched
dvMatrices = cell(1,length(data));
structFieldInds = false(1,length(fields));
for i = 1:length(fields); structFieldInds(i) = isstruct(data(1).(fields{i})); end
structFieldInds = find(structFieldInds);

% get inds of vars in data
varInds = find(ismember({vars.name}, fields));
conditionalInds = find(ismember({conditionals.name}, fields));
dvInds = find(ismember(dvs, fields));

% initialize dvMatrix if dv is found
numDvsFound = sum(isDvFound) + length(dvInds);
if numDvsFound==2
    dvMatrixSize = num2cell([cellfun(@length, {vars.levels}) 2 length(data)]);
    dvMatrix = nan(dvMatrixSize{:});
end
isDvFound(dvInds) = true;



    
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

    % add vars to conditionInds
    skipRow = false;
    for j = varInds
        rowConditionInd = find(ismember(vars(j).levels, data(i).(vars(j).name)));
        if ~isempty(rowConditionInd)
            conditionInds(j) = rowConditionInd;
        else
            skipRow = true;
        end
    end


    if ~skipRow
        
        for j = dvInds; dvsFound(j) = data(i).(dvs{j}); end
        
        if numDvsFound==2
            for j = 1:2
                dvMatrixInds = num2cell([conditionInds j i]);
                dvMatrix(dvMatrixInds{:}) = dvsFound(j);
            end
        
        % loop through nested structs if all dvs have yet to be found
        else
            for j = structFieldInds
                dvMatrices{i} = getSlopeMatrix(data(i).(fields{j}), ...
                    dvs, vars, varsToAvg, varsToSlope, conditionals, conditionInds, dvsFound, isDvFound);
            end

            % avg matrix if data contains field in varsToAvg
            if any(ismember(varsToAvg, fields))
                dvMatrices{i} = nanmean(dvMatrices{i}, dimNum);
            end 
            
            % compute slope
            if any(ismember(varsToSlope, fields))
                
                % iterate through all conditions, replacing second to last dimension with slope of line
                condDims = cellfun(@length, {vars.levels});
                newDvMatrixDims = num2cell([condDims 1 1]);
                newDvMatrix = nan(newDvMatrixDims{:});
                
                for j = 1:prod(cellfun(@length, {vars.levels}))
                    
                    % get inds for regression in dvMatrices{i}
                    condInds = cell(size(condDims));
                    [condInds{:}] = ind2sub(condDims, j);
                    condDimsSub1 = [condInds, {1}, {':'}];
                    condDimsSub2 = [condInds, {2}, {':'}];
                    
                    % get slope relating regressors
                    x = squeeze(dvMatrices{i}(condDimsSub1{:}));
                    y = squeeze(dvMatrices{i}(condDimsSub2{:}));
                    bins = ~isnan(x) & ~isnan(y);
                    fit = polyfit(x(bins), y(bins), 1);
                    
                    % store in newDvMatrix, which has 1 less dim than dvMatrices{i}
                    newDvMatrixInds = [condInds {1} {1}];
                    if sum(bins)>1
                        newDvMatrix(newDvMatrixInds{:}) = fit(1);
                    else
                        newDvMatrix(newDvMatrixInds{:}) = nan;
                    end
                end
                dvMatrices{i} = newDvMatrix;
            end
        end
    end
end


% concatenate matrices obtained from nested structs
if numDvsFound<2; dvMatrix = cat(dimNum, dvMatrices{:}); end

if isFirstPass; dvMatrix = squeeze(dvMatrix); end



