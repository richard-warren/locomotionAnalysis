% function barPlotWrapper(data, dv, vars, varLevels, varLevelNames, varsToAvg)

% to do: get rid of unused fields // collapse across paws, or lower
% levels... // add smp (mouse) name // use all levels of var if not
% explicitly passed in...

% temp
data = flat;
dv = 'isTrialSuccess';
vars = {'isLightOn', 'condition'};
varLevels = {[0,1], {'saline', 'muscimol'}};
varLevelNames = {{'light off', 'light on'}, {'mus', 'sal'}};
varsToAvg = {'trial', 'session', 'mouse'};


% initialiations
if isstruct(data); data = struct2table(data); end
varLevelNums = cellfun(@length, varLevels); % number of levels for each variable
dvMatrix = nan([varLevelNums height(data)]); % each dimension is a variable, and last dimension is data to plotted (eg mouse, trial), which has a length that is determined on the fly
varLevelInds = cell(1,length(vars)); % for every condition (intersection of var levels), this will store the levels of each var
maxConditionSmps = 0;


%% turn all variables into ints, which are easier to work with
fields = data.Properties.VariableNames;
allVarLevels = cell(1,length(fields));
[~,~,varInds] = intersect(vars, fields, 'stable');
[~,~,varsToAvgInds] = intersect(varsToAvg, fields, 'stable');

for i = 1:length(fields)
    if iscell(data.(fields{i})) % if field contains strings    
        [allVarLevels{i}, ~, data.(fields{i})] = unique(data.(fields{i}));
        
        % replace level names with level inds
        varBin = ismember(vars, fields{i});
        if any(varBin)
            [~, ~, varLevels{varBin}] = intersect(varLevels{varBin}, allVarLevels{i}, 'stable');
            varLevels{varBin} = varLevels{varBin}';
        end
    else
        allVarLevels{i} = unique(data.(fields{i}));
    end
end


%%


conditions = combvec(varLevels{:}); % each column is combination of conditions

for i = 1:size(conditions,2)
    
    % get condition bins
    bins = all(table2array(data(:, varInds)) == conditions(:,1)', 2);
    dataSub = data(bins,:);
    
    % avg across varsToAvg
    
    
    
    % save condition data
    indsString = strrep(num2str(cell2mat(varLevelInds)), '  ', ',');
    eval(['dvMatrix(' indsString ',1:length(dataSub)) = [dataSub.(dv)];']);
    if length(dataSub)>maxConditionSmps; maxConditionSmps=sum(dataSub); end
end

eval(['dvMatrix = dvMatrix(' repmat(':,',1,length(vars)) '1:maxConditionSmps);']) % shorten last dimension



% call bar plot function
barPlotRick(dvMatrix, varLevelNames, dv)





