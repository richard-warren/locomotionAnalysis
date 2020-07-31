function [entropies, distributions] = plotEntropies(flat, varargin)

% plot entropy of landing position of first modified paw

% settings
s.condition = '';  % name of field in 'data' that contains different experimental conditions
s.levels = {''};  % levels of s.condition to plot
s.colors = [];
s.binNum = 200;  % number of bins in the distribution
s.poolMice = false;  % whether to pool distributions across all mice

s.xLims = [1 99];  % expressed as percentiles
s.xLimsAbs = [];   % x limits in real world units (m) overwrites xLims if provided

s.successOnly = false;  % whether to only include successful trials
s.modPawOnlySwing = false;  % if true, only include trials where the modified paw is the only one in swing
s.lightOffOnly = false;  % whether to restrict to light on trials
s.deltaMin = 0;  % exclude little step trials where modPawDeltaLength is less than deltaLim standard deviations (standard deviations computed wrt preModPawDeltaLength, which gives a sense of the deltas associated with non-modified steps, which should reflect prediction error more or less...)
s.modSwingContactsMax = 0;  % exclude trials where number of contacts of first mod paw during first mod swing is greater than this value

s.barProps = {};  % properties to pass to barFancy
s.saveLocation = '';  % if provided, save figure automatically to this location

s.verbose = true;


% initialization
if exist('varargin', 'var'); for i = 1:2:length(varargin); s.(varargin{i}) = varargin{i+1}; end; end  % reassign settings passed in varargin
if isempty(s.colors); s.colors = jet(length(s.levels)); end
if isstruct(flat); flat = struct2table(flat); end
cNum = length(s.levels);  % total number of conditions
if s.poolMice; flat.mouse = repmat({''}, height(flat),1); end
mice = unique(flat.mouse);


% restrict to desired trials
if length(s.levels)>1; flat = flat(ismember(flat.(s.condition), s.levels), :); end % keep only trials within condition
if s.successOnly; flat = flat(flat.isTrialSuccess==1, :); end
if s.lightOffOnly; flat = flat(flat.isLightOn==0, :); end
if s.modPawOnlySwing; flat = flat(flat.modPawOnlySwing==1, :); end
if s.modSwingContactsMax
    if s.verbose; fprintf('%.2f of trials removed with modSwingContactsMax criterion\n', ...
            1-mean(flat.modSwingContacts<=s.modSwingContactsMax)); end
    flat = flat(flat.modSwingContacts<=s.modSwingContactsMax, :);
end

if s.deltaMin
    % 1) control distribution, std, all steps
%     minDif = std(flat.preModPawDeltaLength) * 1;
%     bins = abs(flat.modPawDeltaLength)>minDif;
    
    % 2) mod distribution, std, all steps
%     bins = ~(abs(zscore(flat.modPawDeltaLength))<s.deltaMin);
    
    % 3) mod distribution, std, little steps only (old version)
%     bins = ~(abs(zscore(flat.modPawDeltaLength))<s.deltaMin & [flat.isBigStep]==0);

    % 4) mod distribution, meters, little steps only
    bins = ~(abs(flat.modPawDeltaLength)<s.deltaMin & [flat.isBigStep]==0);
    
    % 5) control distribution, std, little steps only per mouse
%     bins = true(1,height(flat));
%     for i = 1:length(mice)
%         mouseBins = strcmp(flat.mouse, mice{i});
%         minDif = std(flat(mouseBins,:).preModPawDeltaLength) * s.deltaMin;
%         if s.verbose; fprintf('%s: minDif %.1f mm\n', mice{i}, minDif*1000); end
%         bins(mouseBins & abs(flat.modPawDeltaLength)<minDif & flat.isBigStep==0) = false;
%     end
    
    if s.verbose; fprintf('%.2f of trials removed with deltaMin criterion\n', 1-mean(bins)); end
    flat = flat(bins,:);
end



% turn the 'condition' into numbers
if ~isempty(s.condition)
    [~, condition] = ismember(flat.(s.condition), s.levels); 
else
    condition = ones(height(flat), 1);
end


% loop over mice
entropies = nan(2, cNum, length(mice));  % (bigstep vs littlestep) X (condition) X (shuffled vs. non-shuffled) X (mice)
distributions = nan(cNum, length(mice), s.binNum);

lims = prctile(flat.modPawDistanceToObs, s.xLims);
lims = lims + [-1 1]*.2*range(lims);
if ~isempty(s.xLimsAbs); lims = s.xLimsAbs; end

x = linspace(lims(1), lims(2), s.binNum);

for h = 1:2
    stepBins = flat.isBigStep==(h-1);
    
    for i = 1:length(mice)
        mouseBins = strcmp(flat.mouse, mice{i});

        % models per condition
        for j = 1:cNum
            conditionBins = condition==j;
            bins = stepBins & mouseBins & conditionBins;
            kd = ksdensity(flat.modPawDistanceToObs(bins), x);
            entropies(h,j,i) = entropy(kd/max(kd));
            
            % compute the overall distribution (pooled across step types)only once
            if h==1
                kd = ksdensity(flat.modPawDistanceToObs(mouseBins & conditionBins), x);
                distributions(j,i,:) = kd;
            end
        end
    end
end



% plot accuracies
figure('position', [200 703.00 300 255.00], 'color', 'white', 'menubar', 'none')
barFancy(entropies, 'ylabel', 'entropy (bits)', 'levelNames', {{'little step', 'big step'}, s.levels}, ...
    'colors', repmat(s.colors,2,1), s.barProps{:})
if ~isempty(s.saveLocation); saveas(gcf, s.saveLocation, 'svg'); end

% plot distributions
figure('position', [206 429 473 120], 'color', 'white', 'menubar', 'none'); hold on
for i = 1:cNum
    plot(x, squeeze(mean(distributions(i,:,:),2)), 'color', s.colors(i,:), 'LineWidth', 2)
end
set(gca, 'YColor', 'none', 'XLim', lims)
plot([0 0], get(gca, 'ylim'), 'color', [0 0 0 .4])
xlabel('landing position (m)')
