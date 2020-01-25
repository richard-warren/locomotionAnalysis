function plotDecisionThresholds(data, varargin)

% settings
s.condition = '';  % name of field in 'data' that contains different experimental conditions
s.levels = {''};  % levels of s.condition to plot
s.colors = [];

s.successOnly = false;  % whether to only include successful trials
s.modPawOnlySwing = false;  % if true, only include trials where the modified paw is the only one in swing
s.lightOffOnly = false;  % whether to restrict to light on trials


% initializations
if exist('varargin', 'var'); for i = 1:2:length(varargin); s.(varargin{i}) = varargin{i+1}; end; end  % reassign settings passed in varargin
if isempty(s.colors); s.colors = jet(length(s.levels)); end

flat = flattenData(data, {'mouse', 'modPawPredictedDistanceToObs', 'modPawDistanceToObs', 'isBigStep', ...
    'isLightOn', 'isTrialSuccess', 'modPawOnlySwing', s.condition});

% restrict to desired trials
if s.successOnly; flat = flat([flat.isTrialSuccess]); end
if s.modPawOnlySwing; flat = flat([flat.modPawOnlySwing]==1); end
if s.lightOffOnly; flat = flat(~[flat.isLightOn]); end

[~, condition] = ismember({flat.(s.condition)}, s.levels);  % turn the 'condition' into numbers
mice = unique({flat.mouse});



thresholds = nan(length(s.levels), length(mice));

for i = 1:length(mice)
    for j = 1:length(s.levels)
        
        bins = strcmp({flat.mouse}, mice{i}) & condition==j;
        
        x = [flat(bins).modPawPredictedDistanceToObs] * 1000;
        y = [flat(bins).isBigStep];
        
        glm = fitglm(x', y', 'Distribution', 'binomial');
        coeffs = glm.Coefficients.Estimate;
        thresholds(j,i) = (-coeffs(1)) / coeffs(2); % solve for prediction = 0
    end
end

figure('position', [2018.00 100.00 252.00 239.00], 'color', 'white', 'menubar', 'none')
barFancy(thresholds, 'ylabel', 'decision threshold', 'levelNames', {s.levels}, 'colors', s.colors)


