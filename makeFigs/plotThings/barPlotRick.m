function barPlotRick(data, conditionNames, dvName, isWithinSubs, smpNames)

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
showScatter = true;
showErrorBars = true;
showViolins = false;
showStats = true;
compareToFirstOnly = length(conditionNames{end})>4; % only run stats between first and all other conditions

groupSeparation = 1;
circSize = 40;
circAlpha = .8;
lineWidth = 1;
lineThickness = 3;
pThresh = .05;

% initializations
condLevels = cellfun(@length, conditionNames);
totalConditions = prod(condLevels);
conditionsMat = nan(length(conditionNames), totalConditions);
labelVertSize = .15*length(conditionNames);
dataDims = size(data);
if length(dataDims)==length(conditionNames); dataDims = [dataDims 1]; end % add singleton dimension if there is only one sample in the dataset
xJitters = linspace(-.5*lineWidth, .5*lineWidth, dataDims(end));
xJitters = xJitters(randperm(length(xJitters)));
if isWithinSubs
    colors = hsv(dataDims(end)) * .8;
else
    colors = hsv(length(conditionNames{end})) * .8;
end

% close all;
% figure('color', 'white', 'menubar', 'none', 'position', [2000 200 interp1([1 16], [200 1200], totalConditions) 300])

% create matrix where each column is an interection of conditions
xPositions = 1:totalConditions;
for i = 1:length(conditionNames)
    repeats = prod(condLevels(i+1:end));
    copies = totalConditions / (repeats*condLevels(i));
    conditionsMat(i,:) = repmat(repelem(1:condLevels(i), repeats), 1, copies);
    xPositions = xPositions + (repelem(1:copies*condLevels(i), repeats)-1) * groupSeparation;
end

% add lines connecting same sample across conditions
if isWithinSubs && dataDims(end)<40 
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
                'linewidth', 1, 'color', [.8 .8 .8]); hold on
        end
    end
end


% plot data
hold on
for i = 1:totalConditions
    inds = cat(1, num2cell(conditionsMat(:,i)), {1:size(data,length(conditionNames)+1)});
    condData = squeeze(data(inds{:}));
    if isWithinSubs; lineColor = [0 0 0]; else; lineColor = colors(conditionsMat(end,i),:); end
    
    % add probability density estimate
    if showViolins
        [p,y] = ksdensity(condData);
        p = p / max(p) * lineWidth*.5; % normalize range
        fill([p -fliplr(p)]+xPositions(i), [y fliplr(y)], [.8 .8 .8], 'FaceColor', 'none')
    end
    
    % scatter raw data
    if showScatter
        if isWithinSubs; scatColor = colors; else; scatColor = [.5 .5 .5]; end
        scatter(xJitters + xPositions(i), condData, ...
            circSize, scatColor, 'filled', 'MarkerFaceAlpha', circAlpha); hold on
    end
    
    if showErrorBars
        err = nanstd(condData);
        line([xPositions(i) xPositions(i)], [err -err] + nanmean(condData), ...
        'color', lineColor, 'linewidth', lineThickness*.5)
    end
    
    % add mean
    line([-.5 .5]*lineWidth + xPositions(i), repmat(nanmean(condData),1,2), ...
        'color', lineColor, 'linewidth', lineThickness)
end
axisData = get(gca);



% add pairwise stats
if showStats && exist('isWithinSubs', 'var')
    [~,~,condInds] = unique(conditionsMat(1:end-1,:)', 'rows'); % only draw lines connecting data across last condition
    
    for i = unique(condInds)'
        dimPairs = nchoosek(1:length(find(condInds==i)), 2); % wrt data dimensions
        condPairs = nchoosek(find(condInds==i), 2); % wrt columns in plot
        
        % sort s.t. more distant comparisons are last
        [~, sortInds] = sort(diff(dimPairs,[],2));
        dimPairs = dimPairs(sortInds,:); condPairs = condPairs(sortInds,:);
        
        % only keep comparisons between first and all other conditions
        if compareToFirstOnly
            bins = any(dimPairs==1,2);
            dimPairs = dimPairs(bins,:);
            condPairs = condPairs(bins,:);
        end
        
        % vertical position of each line, expressed as fraction of y range
        yPosits = linspace(1.0, 1.0+.01*size(dimPairs,1), size(dimPairs,1)); 
        inds = conditionsMat(1:end-1, find(condInds==i,1,'first'))'; % matrix inds for conditions higher up in the hierarchy
        
        for j = 1:size(dimPairs,1)
            
            [inds1, inds2] = deal(cat(2,num2cell([inds dimPairs(j,1)]), {1:size(data,length(conditionNames)+1)}));
            inds2{length(conditionNames)} = dimPairs(j,2);
            
            if isWithinSubs
                [~,p] = ttest(data(inds1{:}), data(inds2{:}));
            else
                [~,p] = ttest2(data(inds1{:}), data(inds2{:}));
            end
            if p<pThresh; lineColor = 'red'; else; lineColor = [.5 .5 .5]; end
            line(xPositions(condPairs(j,:)), max(data(:))*[yPosits(j) yPosits(j)], ...
                'color', lineColor, 'linewidth', 1.0);
        end
    end
    
    set(gca, 'YTick', axisData.YTick, 'YLim', [axisData.YLim(1), max(data(:))*max(yPosits)])
end




% add room beneath x axis for condition labels
yLims = get(gca, 'ylim');
yTicks = get(gca, 'ytick');
set(gca, 'XColor', 'none', ...
    'YLim', [yLims(1)-labelVertSize*range(yLims), yLims(2)], 'YTick', yTicks, ...
    'XLim', [0 xPositions(end)+1], 'XTick', xPositions, 'XTickLabel', [])
line([0 0], [yLims(1)-labelVertSize*range(yLims), yLims(1)], 'color', 'white', 'linewidth', 3) % cover bottom of y axis with white line

% add labels
for i = 1:length(conditionNames)
    
    parentConditions = unique(conditionsMat(1:i-1,:)','rows');
    
    for j = 1:size(parentConditions,1)
        
        if i==1; bins=true(1,totalConditions)'; else; bins = ismember(conditionsMat(1:i-1,:)', parentConditions(j,:), 'rows'); end
        
        for k = 1:condLevels(i)
            inds = find(conditionsMat(i,:)==k & bins');
            xPos = mean(xPositions(inds));
            yPos = yLims(1)-labelVertSize*range(yLims) + ((labelVertSize*range(yLims))/(length(conditionNames)+1))*i;
            if i==length(conditionNames); rotation = 25; else; rotation = 0; end
            condText = text(xPos, yPos, conditionNames{i}(k), 'HorizontalAlignment', 'center', 'rotation', rotation);
            
            % add lines on the side of condition name
            if i<length(conditionNames) && length(conditionNames{i})>1
                textPos = get(condText, 'Extent');
                line([xPositions(inds(1)) textPos(1)], [yPos yPos], 'color', [.5 .5 .5]) % left side of text
                line([textPos(1)+textPos(3) xPositions(inds(end))], [yPos yPos], 'color', [.5 .5 .5]) % right side of text
            end
        end
    end
end
ylabel(dvName)



% add legend
if exist('smpNames', 'var')
    for i = 1:length(smpNames); scatters(i) = scatter(nan,nan,50,colors(i,:),'o','filled'); end % create dummy scatters
    legend(scatters, smpNames, 'Location', 'northeastoutside', 'box', 'off')
end




