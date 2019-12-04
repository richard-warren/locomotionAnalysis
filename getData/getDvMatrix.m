function dvMatrix = getDvMatrix(data, dv, vars, varsToAvg, conditionals, conditionInds)

% given a 'data' struct (see 'getExperimentData'), recursively digs into
% the struct until it finds 'dv', the dependent varialbe of interest //
% returns a matrix for dv broekn down by the varialbes ('vars'), each of
% which has a certain number of 'levels' (e.g. the variable 'sex' has the
% levels 'male', 'female')

% dvMatrix has one dimension per var, and the length of that dimension is
% the number of levels for the var // the final dimenion is the 'subject'
% dimension, and has the dv value for all subjects in a certain condition
% (a condition being an intesection of variable3 levels) // 'conditionals'
% restrict the data to those that satisfy certain conditions

% 'vars' is a struct array with fields 'name' and 'levels'
% 'conditionals' is a struct array with fields 'name' and 'condition'
% 'varsToAvg' is a list of levels in the hierarchy to average over (e.g.
% average the dv within sessions, mice)

% don't use 'conditionInds', which is only used internally by recursive
% calls of the function to keep track of the condition indices of the
% current data... // data at the second level of the first var, and the
% third level of the second var will have indices [2,3] in dvMatrix //
% 'conditionInds' keeps track of this, and is passed down through recursive
% calls so when the dv is found, we know where to put it in dvMatrix


% EXAMPLE
% dvMatrix = (data, ...
%             'trialVel', ...
%             struct('name', 'isLeading', 'levels', [1 0], 'levelNames', {{'leading', 'trailing'}}), ...
%             {'mouse', 'session'}, ...
%             struct('name', 'isLightOn', 'condition', @(x) x==1)

% todo: multiple dvs // how to bin variables // meaningful error messages // make 'mean' or 'median' option


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
    structInds = false(1,length(fields));
    for i = 1:length(fields); structInds(i) = isstruct(data(1).(fields{i})); end
    structInds = find(structInds);
    dvMatrices = cell(1,length(data));
end

% get inds of vars in data
varInds = find(ismember({vars.name}, fields));  % inds of vars that are colums in data
avgDvs = any(ismember(varsToAvg, fields));  % whether we should be averaging the dv at this level of the hierarchy
conditionalInds = find(ismember({conditionals.name}, fields));  % inds of conditionals present at this level of the hierarchy

% return empty matrix if branch contains neither dv nor nested structs
if ~dvFound && isempty(structInds)
    dvMatrix = [];

% otherwise, loop through rows
else
    % find rows that satisfy conditionals
    binsToAnalyze = true(1,length(data));
    if ~isempty(conditionalInds)    
        for i = conditionalInds
            if any(cellfun(@ischar, ({data.(conditionals(i).name)})))  % if contains strings
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
            rowLevelInd = find(ismember(vars(j).levels, data(i).(vars(j).name)));  % ind of level within var for this row (what level of the current variable are we at?)
            if ~isempty(rowLevelInd)
                conditionInds(j) = rowLevelInd;
            else
                skipRow = true;  % skip row if the requested level of any single variable is not present
            end
        end
        
        if ~skipRow
            
            % if dv is in 'data', add to correct location in dvMatrix
            if dvFound
                dvMatrixInds = num2cell([conditionInds i]);
                dvMatrix(dvMatrixInds{:}) = data(i).(dv);

            % otherwise loop through nested structs
            else
                for j = structInds
                    dvMatrices{i} = getDvMatrix(data(i).(fields{j}), ...
                        dv, vars, varsToAvg, conditionals, conditionInds);
                    if ~isempty(dvMatrices{i}); break; end % terminate loop once a branch containing dv is found
                end

                if avgDvs; dvMatrices{i} = nanmean(dvMatrices{i}, length(vars)+1); end  % avg matrix if data contains field in varsToAvg
            end
        end
    end
    
    if ~dvFound
        % todo: when nested structs have different lengths of the final dimension, concatenating in this way will result in many reduntant NaNs...
        dvMatrix = cat(length(vars)+1, dvMatrices{:});  % concatenate matrices obtained from nested structs
    end
end



