function plotDecisionThresholds(flat, varargin)

% settings
s.condition = '';  % name of field in 'data' that contains different experimental conditions
s.levels = {''};  % levels of s.condition to plot
s.colors = [];
s.outcome = 'isBigStep';

s.successOnly = false;  % whether to only include successful trials
s.modPawOnlySwing = false;  % if true, only include trials where the modified paw is the only one in swing
s.lightOffOnly = false;  % whether to restrict to light on trials
s.deltaMin = 0;  % exclude trials where little step length is adjusted less than deltaMin


s.barProps = {};  % properties to pass to barFancy
s.saveLocation = '';  % if provided, save figure automatically to this location


% initializations
if exist('varargin', 'var'); for i = 1:2:length(varargin); s.(varargin{i}) = varargin{i+1}; end; end  % reassign settings passed in varargin
if isempty(s.colors); s.colors = jet(length(s.levels)); end
if isstruct(flat); flat = struct2table(flat); end


% restrict to desired trials
if s.successOnly; flat = flat(flat.isTrialSuccess,:); end
if s.modPawOnlySwing; flat = flat(flat.modPawOnlySwing==1,:); end
if s.lightOffOnly; flat = flat(~flat.isLightOn,:); end
% if s.deltaMin
%     minDif = std(flat.preModPawDeltaLength) * s.deltaMin;
%     flat = flat(abs(flat.modPawDeltaLength)>minDif,:);
% end
if s.deltaMin; flat = flat( ~(abs(zscore(flat.modPawDeltaLength))<s.deltaMin & [flat.isBigStep]==0), :); end

if ~isempty(s.condition)  % if no condition provided, create dummy variable
    [~, condition] = ismember(flat.(s.condition), s.levels);  % turn the 'condition' into numbers
else
    condition = ones(height(flat),1);
end
mice = unique(flat.mouse);


thresholds = nan(length(s.levels), length(mice));

for i = 1:length(mice)
    for j = 1:length(s.levels)
        
        bins = strcmp(flat.mouse, mice{i}) & condition==j;
        
        x = flat.modPawPredictedDistanceToObs(bins) * 1000;
        y = flat.(s.outcome)(bins);
        
        glm = fitglm(x', y', 'Distribution', 'binomial');
        coeffs = glm.Coefficients.Estimate;
        thresholds(j,i) = (-coeffs(1)) / coeffs(2); % solve for prediction = 0
    end
end

figure('position', [2018.00 100.00 252.00 239.00], 'color', 'white', 'menubar', 'none')
barFancy(thresholds, 'ylabel', 'decision threshold', 'levelNames', {s.levels}, 'colors', s.colors, s.barProps{:})


% save
if ~isempty(s.saveLocation); saveas(gcf, s.saveLocation, 'svg'); end


