

% input would be DV, levels of conditions to plot, ad variables to avg over

% !!! first would delete columns not used...

data = struct2table(kinData);

%%


% the following collapses across varsToAvgOver
varsToAvgOver = {'session'};
dv = 'vel';
% conditionBins = false(1,height(data)); conditionBins(conditionInds) = true; % temp
conditionBins = true(1,height(data));

for var = varsToAvgOver
    
    varLevels = unique(data.(var{1}))';
    
    for lev = varLevels
%         keyboard
        indsSub = find(strcmp(data.(var{1}), lev{1}) & conditionBins');
        meanSub = mean(data.(dv)(indsSub));
        data.(dv)(indsSub(1)) = meanSub; % replace first entry with mean
        data(indsSub(2:end),:) = [];
        conditionBins(indsSub(2:end)) = [];
    end 
end