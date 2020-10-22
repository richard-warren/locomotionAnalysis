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
s.ytick = [];
s.addBars = false; % if true, add real bars, instead of just horizontal line at mean
s.barAlpha = .2; % transparency of bars
s.summaryFunction = @nanmean; % user can change this to nanmedian, for example
s.numVariables = [];

s.groupSeparation = 1;
s.circSize = 40;
s.scatAlpha = .6;
s.lineWidth = 1;
s.lineThickness = 3;
s.pThresh = .05;

% reassign settings contained in opts
if exist('opts', 'var'); for i = 1:2:length(opts); s.(opts{i}) = opts{i+1}; end; end

% initializations
figColor = get(gcf, 'color');
if isempty(s.numVariables); s.numVariables = length(size(data))-1; end
varLevelNum = size(data); varLevelNum = varLevelNum(1:s.numVariables); % number of levels for each variable
totalConditions = prod(varLevelNum);
dataDims = size(data);
if s.numVariables==length(dataDims); dataDims = [dataDims 1]; end  % add singleton dimension if there is only one sample per condition

if ischar(s.conditionColors) % set bar colors if color is specified as a string
    s.conditionColors = eval([s.conditionColors '(totalConditions)']);
elseif isequal(size(s.conditionColors), [1 3]) % if specified as a single rbg value, replicate into a matrix
    s.conditionColors = repmat(s.conditionColors,totalConditions,1);
end

if ischar(s.scatColors) % set bar colors if color is specified as a string
    if s.isWithinSubs
        s.scatColors = eval([s.scatColors '(dataDims(end))']);
    else
        s.scatColors = [.5 .5 .5];
    end
end

conditionsMat = nan(s.numVariables, totalConditions);
labelVertSize = .15*s.numVariables; % size of space below figure to give to to axis labels, expressed as fraction of y range
statsVertSpacing = .02; % vertical spacing of stat comparison lines, expressed as fraction of y range
xJitters = linspace(-.5*s.lineWidth, .5*s.lineWidth, dataDims(end));
xJitters = xJitters(randperm(length(xJitters)));
hold on


% create matrix where each column is an interection of conditions
xPositions = 1:totalConditions;
for i = 1:s.numVariables
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
            try
            inds = num2cell([conditionsMat(:,j); i]);
            smpData(j) = squeeze(data(inds{:}));
            catch; keyboard; end
        end
        
        % draw lines connecting data only across levels of last condition
        for j = unique(condInds)'
            line(xPositions(condInds==j)+xJitters(i), smpData(condInds==j), ...
                'linewidth', 1, 'color', [.8 .8 .8 s.scatAlpha]);
        end
    end
end


% plot data
allData = cell(1,totalConditions);
for i = 1:totalConditions
    inds = cat(1, num2cell(conditionsMat(:,i)), {1:size(data,length(dataDims))});
    condData = squeeze(data(inds{:}));
    allData{i} = condData;
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
        line([xPositions(i) xPositions(i)], [err -err] + s.summaryFunction(condData), ...
            'color', s.conditionColors(i,:), 'linewidth', s.lineThickness*.5)
    end
    
    % add mean
    if ~s.addBars
        line([-.5 .5]*s.lineWidth + xPositions(i), repmat(s.summaryFunction(condData),1,2), ...
            'color', s.conditionColors(i,:), 'linewidth', s.lineThickness)
    end
end


% SET Y LIMITS
if isempty(s.ylim)
    s.ylim = get(gca, 'YLim');
    
    if s.addBars  % bar plot will always start at zero // otherwise, use the automatically determined y limits
        
        % collect all data included in plot to determine range
        summaries = cellfun(s.summaryFunction, allData)';
        errors = cellfun(@nanstd, allData)';
        ys = summaries;  % ys contains all data to be included in range, which depends on elements are to be included in plot
        if s.showScatter; ys = [ys; data(:)]; end
        if s.showErrorBars; ys = [ys; summaries+errors; summaries-errors]; end
        
        if min(ys)>0  % if all data are positive, lower y limit is 0
            s.ylim(1) = 0;
        elseif max(ys)<0  % if all data are negative, upper y limit is zero
            s.ylim(2) = 0;
        end
    end
end
set(gca, 'YLim', s.ylim);

line([0 xPositions(end)+s.lineWidth/2], [0 0], 'color', get(gca, 'YColor'))  % add line at y=0 zero

if ~isempty(s.ytick)
    set(gca, 'YTick', s.ytick);
    yTicks = s.ytick;
else
    yTicks = get(gca, 'ytick');
end

yMin = s.ylim(1)-labelVertSize*range(s.ylim);  % this is where new bottom of figure will be after adding room below for labels

% add bars
if s.addBars
    for i = 1:totalConditions
        
        xs = [-.5 .5]*s.lineWidth + xPositions(i);
        ys = [0, s.summaryFunction(allData{i})];

        plot([xs(1) xs(1) xs(2) xs(2)], [ys(1) ys(2) ys(2) ys(1)], ...
            'LineWidth', s.lineThickness, 'Color', s.conditionColors(i,:));
        
        if s.barAlpha>0
            pshape= polyshape([xs(1) xs(1) xs(2) xs(2) xs(1)], [ys(1) ys(2) ys(2) ys(1) ys(1)]);
            pshape = plot(pshape);
            set(pshape, 'EdgeColor', 'none', 'FaceColor', [s.conditionColors(i,:) s.barAlpha])
        end
    end
end



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
        yMax = s.ylim(2) + statsVertSpacing*range(s.ylim)*size(dimPairs,1);
        yPosits = linspace(s.ylim(2), yMax, size(dimPairs,1)); 
        inds = conditionsMat(1:end-1, find(condInds==i,1,'first'))'; % matrix inds for conditions higher up in the hierarchy
        
        for j = 1:size(dimPairs,1)
            
            [inds1, inds2] = deal(cat(2,num2cell([inds dimPairs(j,1)]), {1:size(data,length(dataDims))}));
            inds2{s.numVariables} = dimPairs(j,2);
            
            if s.isWithinSubs
                [~,p] = ttest(data(inds1{:}), data(inds2{:}));
            else
                [~,p] = ttest2(data(inds1{:}), data(inds2{:}));
            end
            if p<s.pThresh; lineColor = 'red'; else; lineColor = [.5 .5 .5]; end
            line(xPositions(condPairs(j,:)), [yPosits(j) yPosits(j)], ...
                'color', lineColor, 'linewidth', 1.0);
        end
    end
end




% cover space above and below s.ylim with boxes to occlude bars outside of limits
xs = [0 xPositions(end)+1];
ys = {[s.ylim(1)-range(s.ylim), s.ylim(1)], ... % bottom box
      [s.ylim(2), s.ylim(2)+range(s.ylim)]};    % top box
for i = ys
%     keyboard
    rectangle('Position', [xs(1) i{1}(1) range(xs) range(i{1})], ...
        'FaceColor', figColor, 'EdgeColor', 'none');
end



% add labels
for i = 1:length(s.conditionNames)
    
    parentConditions = unique(conditionsMat(1:i-1,:)','rows');
    
    for j = 1:size(parentConditions,1)
        
        if i==1; bins=true(1,totalConditions)'; else; bins = ismember(conditionsMat(1:i-1,:)', parentConditions(j,:), 'rows'); end
        
        for k = 1:varLevelNum(i)
            inds = find(conditionsMat(i,:)==k & bins');
            xPos = mean(xPositions(inds));
            yPos = s.ylim(1)-labelVertSize*range(s.ylim) + ((labelVertSize*range(s.ylim))/length(dataDims)*i);
            if i==s.numVariables; rotation = 25; else; rotation = 0; end
            if ~isempty(s.conditionNames)
                condText = text(xPos, yPos, s.conditionNames{i}(k), 'rotation', rotation, ...
                    'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
            end
            
            % add lines on the side of condition name
            if i<s.numVariables
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


% add room above figure for stat lines
if s.showStats
    line([0 0], [s.ylim(2), yMax], 'color', figColor, 'linewidth', 3) % cover top of y axis with white line
    s.ylim = [s.ylim(1), yMax];
end


% add room below figure for labels
line([0 0], [yMin, s.ylim(1)], 'color', figColor, 'linewidth', 3) % cover bottom of y axis with white line
s.ylim = [yMin, s.ylim(2)];

set(gca, 'YLim', s.ylim, 'YTick', yTicks, ...
    'XLim', [0 xPositions(end)+1], 'XColor', 'none', ...
    'Color', figColor)


% add legend
if ~isempty(s.smpNames)
    for i = 1:length(s.smpNames); scatters(i) = scatter(nan,nan,50,s.scatColors(i,:),'o','filled'); end % create dummy scatters
    legend(scatters, s.smpNames, 'Location', 'northeastoutside', 'box', 'off')
end


% add y axis label
if ~isempty(s.ylabel)
    lab = ylabel(s.ylabel);
    labPos = get(lab, 'position');
    labPos(2) = mean(s.ylim);
    set(lab, 'position', labPos);
end

pause(.001) % when doing many subplots, this makes sure they pop up one by one




