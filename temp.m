%% play around with generating propensity scores for motor cortex data

load(fullfile(getenv('OBSDATADIR'), 'matlabData', 'mtc_muscimol_data.mat'), 'data');

%% get full data table
locoVars = {'velAtWiskContact', 'angleAtWiskContactContra', 'tailHgt'};

flat = flattenData(data, [{'mouse', 'session', 'trial', 'condition', 'isTrialSuccess', 'isLightOn'} locoVars]);
flat = struct2table(flat);
varBins = ismember(flat.Properties.VariableNames, locoVars);
[~,~,y] = unique(flat.condition); % condition encoded numerically
isManip = ~logical(y-1);
flat.isManip = isManip;


%% get sub data table

% flatSub = flat(strcmp(flat.mouse, 'sen4') & ~flat.isLightOn, :);
flatSub = flat;

validBins = all(~isnan([table2array(flatSub(:, varBins)), flatSub.isManip]),2);
validBins = ones(1, height(flatSub));
% flatSub = flatSub(validBins,:);

X = table2array(flatSub(:,varBins));
y = flatSub.isManip;


pairs = propensityMatching(X, y, ...
    {'percentileThresh', 100, 'predictorNames', locoVars});

%%
flatMatched = flatSub(pairs(:),:);
fprintf('\nvel decrease: %.2f, success decrease: %.2f\n', ...
    nanmean(flatSub.velAtWiskContact(flatSub.isManip)) - nanmean(flatSub.velAtWiskContact(~flatSub.isManip)), ...
    mean(flatSub.isTrialSuccess(flatSub.isManip)) - mean(flatSub.isTrialSuccess(~flatSub.isManip)));
fprintf('matched vel decrease: %.2f, matched success decrease: %.2f\n', ...
    nanmean(flatMatched.velAtWiskContact(flatMatched.isManip)) - nanmean(flatMatched.velAtWiskContact(~flatMatched.isManip)), ...
    mean(flatMatched.isTrialSuccess(flatMatched.isManip)) - mean(flatMatched.isTrialSuccess(~flatMatched.isManip)));


