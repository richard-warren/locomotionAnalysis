function plotPredictors(flat, predictors, target, varargin)

% plots probability of 'target' being true as a function each 'predictor'
% // plots one row per condition // can pool across mice or average across
% mouse averages



% settings
s.condition = '';  % name of field in 'data' that contains different experimental conditions
s.levels = {''};  % levels of s.condition to plot
s.colors = [];
s.names = {};  % use these names to label subplots (overwriting names in 'predictors')
s.overlayConditions = false;  % whether to plot different conditions on top of one another, or as different rows
s.figPos = [];
s.fontSize = 10;

s.avgMice = false;
s.percentileLims = [10 99];
s.binNum = 100;
s.binWidth = .15;  % expressed in terms of frace of x limits
s.yLim = [0 1];  % if using this for non-binary dependent variable this may become useful
s.mouseAlpha = .25;
s.subplotDims = [];

s.successOnly = false;  % whether to only include successful trials
s.modPawOnlySwing = false;  % if true, only include trials where the modified paw is the only one in swing
s.lightOffOnly = false;  % whether to restrict to light on trials
s.deltaMin = 0;  % exclude little step trials where modPawDeltaLength is less than deltaLim standard deviations
s.modSwingContactsMax = 0;  % exclude trials where number of contacts of first mod paw during first mod swing is greater than this value

s.lineWidth = .05;  % expressed as fraction of y axis
s.mouseLineWidth = .02;  % expressed as fraction of y axis
s.verbose = true;

s.textOnly = false;

s.saveLocation = '';  % if provided, save figure automatically to this location


% initialization
if exist('varargin', 'var'); for i = 1:2:length(varargin); s.(varargin{i}) = varargin{i+1}; end; end  % reassign settings passed in varargin
if length(s.levels)==1; colorNum=length(predictors); else; colorNum=length(s.condition); end 
if isempty(s.colors); s.colors = jet(colorNum); end
if isstruct(flat); flat = struct2table(flat); end
cNum = length(s.levels);  % total number of conditions
if isempty(s.names); s.names = predictors; end
if isempty(s.subplotDims)
    if s.overlayConditions
        s.subplotDims = [1 length(predictors)];
    else
        s.subplotDims = [cNum length(predictors)];
    end
end
if isempty(s.figPos); s.figPos = [100 400 200*s.subplotDims(2) 150*s.subplotDims(1)]; end

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
end
if s.verbose; fprintf('%.2f of trials removed with deltaMin criterion\n', 1-mean(bins)); end
flat = flat(bins,:);


if ~s.avgMice; flat.mouse = repmat({'temp'}, height(flat), 1); end  % if pooling across mice, rename all mice with dummy string
mice = unique(flat.mouse);
if ~isempty(s.condition)
    [~, condition] = ismember(flat.(s.condition), s.levels);  % turn the 'condition' into numbers
else
    condition = ones(height(flat), 1);
end

figure('position', s.figPos, 'color', 'white', 'menubar', 'none')

% loop over conditions
[probs, thicknesses] = deal(zeros(cNum, length(predictors), length(mice), s.binNum));
% 
for j = 1:length(predictors)
    xLims = prctile(flat.(predictors{j}), s.percentileLims);  % base x limits on distribution for this predictor across conditions and across mice
    binCenters = linspace(xLims(1), xLims(2), s.binNum);
    binWidth = range(xLims) * s.binWidth;
    
    for i = 1:cNum
        if s.overlayConditions
            sInd = j;
        else
            sInd = (i-1)*length(predictors) + j;
        end
        
        subplot(s.subplotDims(1), s.subplotDims(2), sInd); hold on
        if length(s.levels)==1; c = s.colors(j,:); else; c=s.colors(i,:); end
        
        for k = 1:length(mice)
            mouseBins = strcmp(flat.mouse, mice{k});
            
            x = flat.(predictors{j})(mouseBins & condition==i);
            y = flat.(target)(mouseBins & condition==i);
            
            for m = 1:s.binNum
                bins = x>(binCenters(m)-binWidth*.5) & x<=(binCenters(m)+binWidth*.5);
                if sum(bins)>0
                    probs(i,j,k,m) = nanmean(y(bins));
                    thicknesses(i,j,k,m) = sum(bins) / length(x);
                end
            end
            
            % plot mouse
            if s.avgMice
                if s.mouseAlpha
                    p = squeeze(probs(i,j,k,:));
                    t = squeeze(thicknesses(i,j,k,:));
                    t = t*s.mouseLineWidth/max(t);
                    if ~s.textOnly; patch([binCenters fliplr(binCenters)]', [p+t; flipud(p-t)], c, 'EdgeColor', 'none', 'FaceAlpha', s.mouseAlpha); end
                end
            end
        end
        
        % plot avg across mice
        p = squeeze(nanmean(probs(i,j,:,:),3));
        t = squeeze(nanmean(thicknesses(i,j,:,:),3));
        t = t*s.lineWidth/max(t);
        if ~s.textOnly; patch([binCenters fliplr(binCenters)], [p+t; flipud(p-t)], c, 'EdgeColor', 'none'); end
        
        if j==1; ylabel(s.levels{i}); end
        if i==cNum; xlabel(s.names{j}); end
        if s.textOnly
            set(gca, 'yLim', s.yLim, 'XTick', [], 'YTick', [0 1])
        else
            set(gca, 'yLim', s.yLim, 'XTick', [], 'xlim', xLims, 'YTick', [0 1], ...
                'FontName', 'Arial', 'FontSize', s.fontSize)
        end
        if mod(j-1, s.subplotDims(2))~=0; set(gca, 'YTickLabel', []);  % only label y axes for left-most plots
    end
end


% save
if ~isempty(s.saveLocation); saveas(gcf, s.saveLocation, 'svg'); end

end
