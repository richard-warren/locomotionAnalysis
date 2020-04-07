function plotModelPredictors(flat, predictors, target, varargin)

% plots probability of 'target' being true as a function each 'predictor'
% // plots one row per condition // can pool across mice or average across
% mouse averages



% settings
s.condition = '';  % name of field in 'data' that contains different experimental conditions
s.levels = {''};  % levels of s.condition to plot
s.colors = [];

s.avgMice = false;
s.percentileLims = [1 99];
s.binNum = 100;
s.binWidth = .1;  % expressed in terms of frace of x limits
s.yLim = [0 1];  % if using this for non-binary dependent variable this may become useful

s.successOnly = false;  % whether to only include successful trials
s.modPawOnlySwing = false;  % if true, only include trials where the modified paw is the only one in swing
s.lightOffOnly = false;  % whether to restrict to light on trials
s.deltaMin = 0;  % exclude little step trials where modPawDeltaLength is less than deltaLim standard deviations

s.saveLocation = '';  % if provided, save figure automatically to this location


% initialization
if exist('varargin', 'var'); for i = 1:2:length(varargin); s.(varargin{i}) = varargin{i+1}; end; end  % reassign settings passed in varargin
if length(s.levels)==1; colorNum=length(predictors); else; colorNum=length(s.condition); end 
if isempty(s.colors); s.colors = jet(colorNum); end
if isstruct(flat); flat = struct2table(flat); end
cNum = length(s.levels);  % total number of conditions

% restrict to desired trials
if s.successOnly; flat = flat(flat.isTrialSuccess==1, :); end
if s.lightOffOnly; flat = flat(flat.isLightOn==0, :); end
if s.modPawOnlySwing; flat = flat(flat.modPawOnlySwing==1, :); end
if s.deltaMin
    minDif = std(flat.preModPawDeltaLength) * s.deltaMin;
    flat = flat(abs(flat.modPawDeltaLength)>minDif,:);
end

if ~s.avgMice; flat.mouse = repmat({'temp'}, height(flat), 1); end  % if pooling across mice, rename all mice with dummy string
mice = unique(flat.mouse);
if ~isempty(s.condition)
    [~, condition] = ismember(flat.(s.condition), s.levels);  % turn the 'condition' into numbers
else
    condition = ones(height(flat), 1);
end

figure('position', [100 400 200*length(predictors) 150*length(s.levels)], ...
    'color', 'white', 'menubar', 'none')

% loop over conditions
[probs, transparencies] = deal(nan(cNum, length(predictors), length(mice), s.binNum));

for j = 1:length(predictors)
    xLims = prctile(flat.(predictors{j}), s.percentileLims);  % base x limits on distribution for this predictor across conditions and across mice
    binCenters = linspace(xLims(1), xLims(2), s.binNum);
    binWidth = range(xLims) * s.binWidth;
    
    for i = 1:cNum
        
        subplot(cNum, length(predictors), (i-1)*length(predictors) + j); hold on
        if length(s.levels)==1; c = s.colors(j,:); else; c=s.colors(i,:); end
        
        for k = 1:length(mice)
            mouseBins = strcmp(flat.mouse, mice{k});
            
            x = flat.(predictors{j})(mouseBins & condition==i);
            y = flat.(target)(mouseBins & condition==i);
            
            for m = 1:s.binNum
                bins = x>binCenters(m)-binWidth*.5 & x<=binCenters(m)+binWidth*.5;
                probs(i,j,k,m) = nanmean(y(bins));
                transparencies(i,j,k,m) = sum(bins) / length(y);
            end
            
            % plot mouse
            if s.avgMice
                x = squeeze(probs(i,j,k,:));
                t = squeeze(transparencies(i,j,k,:));
                patch(binCenters, [x(1:end-1)' NaN], c, 'EdgeColor', c, ...
                    'FaceVertexAlphaData', t/max(t)*.5, ...
                    'AlphaDataMapping', 'none', 'EdgeAlpha', 'interp', 'LineWidth', 1)
            end
        end
        
        % plot avg across mice
        x = squeeze(nanmean(probs(i,j,:,:),3));
        t = squeeze(nanmean(transparencies(i,j,:,:),3));
        patch(binCenters, [x(1:end-1)' NaN], c, 'EdgeColor', c, ...
            'FaceVertexAlphaData', t/max(t), ...
            'AlphaDataMapping', 'none', 'EdgeAlpha', 'interp', 'LineWidth', 4)
        
        if j==1; ylabel(s.levels{i}); end
        if i==1; title(predictors{j}); end
        set(gca, 'yLim', s.yLim, 'XTick', [])
    end
end



% save
if ~isempty(s.saveLocation); saveas(gcf, s.saveLocation, 'svg'); end

end
