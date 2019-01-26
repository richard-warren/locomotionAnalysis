function barPlotRick(data, conditionNames, dvName, smpNames)

% TO DO: add ability to make arbitrary statistical comparisons, within and
% between subs // add box plots and/or violin plot options // document //
% get rid of eval statements // replace my hacky code with combvec!

% temp
% conditionNames = {{'light on', 'light off'}, {'foreLead', 'foreLag', 'hindLead', 'hindLab'}, {'ipsi', 'contra'}, {'sal', 'mus'}};
% smpNames = {'run1', 'run2', 'run3', 'run4', 'run5'};
% mouseNum = 5;
% dv = 'vel';
% data = rand([cellfun(@length, conditionNames), mouseNum]);
% data(:,:,2,:) = data(:,:,2,:) * .5;
% data(:,2,:,:) = data(:,2,:,:) * .75;

% settings
connectLines = false;
groupSeparation = 1.5;
circSize = 40;
circAlpha = .8;
lineWidth = 1;

% initializations
condLevels = cellfun(@length, conditionNames);
totalConditions = prod(condLevels);
conditionsMat = nan(length(conditionNames), totalConditions);
labelVertSize = .1*length(conditionNames);
dataDims = size(data);
colors = hsv(dataDims(end)) * .8;
xJitters = linspace(-.5*lineWidth, .5*lineWidth, dataDims(end));
xJitters = xJitters(randperm(length(xJitters)));
conditionSymbols = {'o', 's', 'h', 'd'};

close all;
figure('color', 'white', 'menubar', 'none', 'position', [0 194 1345 314])

% create matrix where each column is an interection of conditions
xPositions = 1:totalConditions;
for i = 1:length(conditionNames)
    repeats = prod(condLevels(i+1:end));
    copies = totalConditions / (repeats*condLevels(i));
    conditionsMat(i,:) = repmat(repelem(1:condLevels(i), repeats), 1, copies);
    xPositions = xPositions + (repelem(1:copies*condLevels(i), repeats)-1) * groupSeparation;
end

% add lines connecting same sample across conditions
if connectLines
    for i = 1:dataDims(end)
        smpData = nan(1,totalConditions);
        for j = 1:totalConditions
            indsString = strrep(num2str(conditionsMat(:,j)'), '  ', ',');
            smpData(j) = squeeze(eval(['data(' indsString ',i)']));
        end
        line(xPositions+xJitters(i), smpData, ...
            'linewidth', 1, 'color', [.8 .8 .8]); hold on
    end
end


% plot data
hold on
for i = 1:totalConditions
    indsString = strrep(num2str(conditionsMat(:,i)'), '  ', ',');
    condData = squeeze(eval(['data(' indsString ',:)']));
    
    % add probability density estimate
    [p,y] = ksdensity(condData);
    p = p / max(p) * lineWidth*.5;
    fill([p -fliplr(p)]+xPositions(i), [y fliplr(y)], [.8 .8 .8], 'FaceColor', 'none')
    
    % scatter raw data
    scatter(xJitters + xPositions(i), condData, ...
        circSize, colors, conditionSymbols{conditionsMat(end,i)}, 'filled', 'MarkerFaceAlpha', circAlpha); hold on
    
    % add mean
    line([-.5 .5]*lineWidth + xPositions(i), repmat(mean(condData),1,2), ...
        'color', 'black', 'linewidth', 2)
    
    
end

% add room beneath x axis for condition labels
yLims = get(gca, 'YLim');
yTicks = get(gca, 'YTick');
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
            condText = text(xPos, yPos, conditionNames{i}(k), 'HorizontalAlignment', 'center');
            
            % add lines on the side of condition name
            if i<length(conditionNames)
                textPos = get(condText, 'Extent');
                line([xPositions(inds(1)) textPos(1)], [yPos yPos], 'color', [.5 .5 .5]) % left side of text
                line([textPos(1)+textPos(3) xPositions(inds(end))], [yPos yPos], 'color', [.5 .5 .5]) % right side of text
            end
        end
    end
end
ylabel(dvName, 'position', [-.8 mean(yLims) -1])

% add legend
if exist('smpNames', 'var')
    for i = 1:length(smpNames); scatters(i) = scatter(nan,nan,50,colors(i,:),'o','filled'); end % create dummy scatters
    legend(scatters, smpNames, 'Location', 'northeastoutside', 'box', 'off')
end




