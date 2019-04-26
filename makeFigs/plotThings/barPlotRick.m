function barPlotRick(data, opts)

% TO DO: should allow user to define whether each var is between or within, and also to define arbitrary stat comparisons //deal with legend...
% document // replace my hacky code with combvec! //

% % create fake data
% conditionNames = {{'light on', 'light off'}, {'foreLead', 'foreLag', 'hindLead', 'hindLab'}, {'ipsi', 'contra'}, {'sal', 'mus'}};
% smpNames = {'run1', 'run2', 'run3', 'run4', 'run5'};
% mouseNum = 5;
% dv = 'vel';
% data = rand([cellfun(@length, conditionNames), mouseNum]);
% data(:,:,2,:) = data(:,:,2,:) * .5;
% data(:,2,:,:) = data(:,2,:,:) * .75;

% settings
s.showScatter = true;
s.showErrorBars = true;
s.showViolins = false;
s.showStats = true;
s.compareToFirstOnly = true; % only run stats between first and all other conditions
s.conditionNames = {}; % nested cell array, where each cell array contains names of the levels of each variable
s.isWithinSubs = true; % are the conditions within subjects
s.smpNames = {};
s.ylabel = [];
s.conditionColors = [.2 .2 .2];
s.violinAlpha = .2; % alpha for violin plot fill
s.scatColors = 'hsv'; % can be a single rgb value, a name of a color space, or a matrix of colors
s.ylim = [];

s.groupSeparation = 1;
s.circSize = 40;
s.scatAlpha = .6;
s.lineWidth = 1;
s.lineThickness = 3;
s.pThresh = .05;

% reassign settings contained in opts
if exist('opts', 'var'); for i = 1:2:length(opts); s.(opts{i}) = opts{i+1}; end; end

% initializations
numVariables = length(size(data))-1;
varLevelNum = size(data); varLevelNum = varLevelNum(1:end-1); % number of levels for each variable
totalConditions = prod(varLevelNum);
dataDims = size(data);

if ischar(s.conditionColors) % set bar colors if color is specified as a string
    s.conditionColors = eval([s.conditionColors '(totalConditions)']);
elseif isequal(size(s.conditionColors), [1 3]) % if specified as a single rbg value, replicate into a matrix
    s.conditionColors = repmat(s.conditionColors,totalConditions,1);
end

if ischar(s.scatColors) % set bar colors if color is specified as a string
    if s.isWithinSub
        s.scatColors = eval([s.scatColors '(dataDims(end))']);
    else
        s.scatColors = [.5 .5 .5];
    end
end

conditionsMat = nan(numVariables, totalConditions);
labelVertSize = .15*numVariables;
xJitters = linspace(-.5*s.lineWidth, .5*s.lineWidth, dataDims(end));
xJitters = xJitters(randperm(length(xJitters)));
hold on


% create matrix where each column is an interection of conditions
xPositions = 1:totalConditions;
for i = 1:numVariables
    repeats = prod(varLevelNum(i+1:end));
    copies = totalConditions / (repeats*varLevelNum(i));
    conditionsMat(i,:) = repmat(repelem(1:varLevelNum(i), repeats), 1, copies);
    xPositions = xPositions + (repelem(1:copies*varLevelNum(i), repeats)-1) * s.groupSeparation;
end

% add lines connecting same sample across conditions
if s.isWithinSubs && dataDims(end)<40 
    [~,~,condInds] = unique(conditionsMat(1:end-1,:)', 'rows'); % only draw lines connecting data across last condition
    
    for i = 1:dataDims(end) % for each subject        
        
        % get data in all conditions
        smpData = nan(1,totalConditions);
        for j = 1:totalConditions
            inds = num2cell([conditionsMat(:,j); i]);
            smpData(j) = squeeze(data(inds{:}));
        end
        
        % draw lines connecting data only across levels of last condition
        for j = unique(condInds)'
            line(xPositions(condInds==j)+xJitters(i), smpData(condInds==j), ...
                'linewidth', 1, 'color', [.8 .8 .8]);
        end
    end
end


% plot data
for i = 1:totalConditions
    inds = cat(1, num2cell(conditionsMat(:,i)), {1:size(data,length(dataDims))});
    condData = squeeze(data(inds{:}));
    if s.isWithinSubs; lineColor = [0 0 0]; else; lineColor = scatColors(conditionsMat(end,i),:); end
    
    % add probability density estimate
    if s.showViolins
        [p,y] = ksdensity(condData);
        p = p / max(p) * s.lineWidth*.5; % normalize range
        fill([p -fliplr(p)]+xPositions(i), [y fliplr(y)], [.8 .8 .8], ...
            'FaceColor', s.conditionColors(i,:), 'FaceAlpha', s.violinAlpha, 'EdgeColor', s.conditionColors(i,:))
    end
    
    % scatter raw data
    if s.showScatter
        scatter(xJitters + xPositions(i), condData, ...
            s.circSize, s.scatColors, 'filled', 'MarkerFaceAlpha', s.scatAlpha); hold on
    end
    
    if s.showErrorBars
        err = nanstd(condData);
        line([xPositions(i) xPositions(i)], [err -err] + nanmean(condData), ...
        'color', s.conditionColors(i,:), 'linewidth', s.lineThickness*.5)
    end
    
    % add mean
    line([-.5 .5]*s.lineWidth + xPositions(i), repmat(nanmean(condData),1,2), ...
        'color', s.conditionColors(i,:), 'linewidth', s.lineThickness)
end
axisData = get(gca);



% add pairwise stats
if s.showStats
    [~,~,condInds] = unique(conditionsMat(1:end-1,:)', 'rows'); % only draw lines connecting data across last condition
    
    for i = unique(condInds)'
        dimPairs = nchoosek(1:length(find(condInds==i)), 2); % wrt data dimensions
        condPairs = nchoosek(find(condInds==i), 2); % wrt columns in plot
        
        % sort s.t. more distant comparisons are last
        [~, sortInds] = sort(diff(dimPairs,[],2));
        dimPairs = dimPairs(sortInds,:); condPairs = condPairs(sortInds,:);
        
        % only keep comparisons between first and all other conditions
        if s.compareToFirstOnly
            bins = any(dimPairs==1,2);
            dimPairs = dimPairs(bins,:);
            condPairs = condPairs(bins,:);
        end
        
        % vertical position of each line, expressed as fraction of y range
        yPosits = linspace(1.0, 1.0+.01*size(dimPairs,1), size(dimPairs,1)); 
        inds = conditionsMat(1:end-1, find(condInds==i,1,'first'))'; % matrix inds for conditions higher up in the hierarchy
        
        for j = 1:size(dimPairs,1)
            
            [inds1, inds2] = deal(cat(2,num2cell([inds dimPairs(j,1)]), {1:size(data,length(dataDims))}));
            inds2{numVariables} = dimPairs(j,2);
            
            if s.isWithinSubs
                [~,p] = ttest(data(inds1{:}), data(inds2{:}));
            else
                [~,p] = ttest2(data(inds1{:}), data(inds2{:}));
            end
            if p<s.pThresh; lineColor = 'red'; else; lineColor = [.5 .5 .5]; end
            line(xPositions(condPairs(j,:)), max(data(:))*[yPosits(j) yPosits(j)], ...
                'color', lineColor, 'linewidth', 1.0);
        end
    end
    
    set(gca, 'YTick', axisData.YTick, 'YLim', [axisData.YLim(1), max(data(:))*max(yPosits)])
end




% add room beneath x axis for condition labels
if ~isempty(s.ylim); set(gca, 'YLim', s.ylim); end
yLims = get(gca, 'ylim');
yTicks = get(gca, 'ytick');
set(gca, 'XColor', 'none', ...
    'YLim', [yLims(1)-labelVertSize*range(yLims), yLims(2)], 'YTick', yTicks, ...
    'XLim', [0 xPositions(end)+1], 'XTick', xPositions, 'XTickLabel', [])
line([0 0], [yLims(1)-labelVertSize*range(yLims), yLims(1)], 'color', 'white', 'linewidth', 3) % cover bottom of y axis with white line

% add labels
for i = 1:numVariables
    
    parentConditions = unique(conditionsMat(1:i-1,:)','rows');
    
    for j = 1:size(parentConditions,1)
        
        if i==1; bins=true(1,totalConditions)'; else; bins = ismember(conditionsMat(1:i-1,:)', parentConditions(j,:), 'rows'); end
        
        for k = 1:varLevelNum(i)
            inds = find(conditionsMat(i,:)==k & bins');
            xPos = mean(xPositions(inds));
            yPos = yLims(1)-labelVertSize*range(yLims) + ((labelVertSize*range(yLims))/length(dataDims)*i);
            if i==numVariables; rotation = 25; else; rotation = 0; end
            if ~isempty(s.conditionNames)
                condText = text(xPos, yPos, s.conditionNames{i}(k), 'HorizontalAlignment', 'center', 'rotation', rotation);
            end
            
            % add lines on the side of condition name
            if i<numVariables
                if ~isempty(s.conditionNames)
                    textPos = get(condText, 'Extent');
                    line([xPositions(inds(1)) textPos(1)], [yPos yPos], 'color', [.5 .5 .5]) % left side of text
                    line([textPos(1)+textPos(3) xPositions(inds(end))], [yPos yPos], 'color', [.5 .5 .5]) % right side of text
                else
                    line([xPositions(inds(1)) xPositions(inds(end))], [yPos yPos], 'color', [.5 .5 .5]) % line spanning variable
                end
            end
        end
    end
end
if ~isempty(s.ylabel); ylabel(s.ylabel); end



% add legend
if ~isempty(s.smpNames)
    for i = 1:length(s.smpNames); scatters(i) = scatter(nan,nan,50,scatColors(i,:),'o','filled'); end % create dummy scatters
    legend(scatters, s.smpNames, 'Location', 'northeastoutside', 'box', 'off')
end




