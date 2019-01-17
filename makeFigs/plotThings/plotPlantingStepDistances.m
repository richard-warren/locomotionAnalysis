% function plotPlantingStepDistances(data, figTitle)

% plots characteristics of distace of step preceding step over obstacle //
% idea is to test whether planting distance to obs is something that is
% actively controlled

% temp
figTitle = 'test';
data = kinData; clear kinData
data = data(~[data.isLightOn]);

% settings
xLims = [-.07 .01];
yLims = xLims;
pawsToPlot = [1 2]; % {'leading fore', 'lagging fore', 'leading hind', 'lagging hind'};
binSize = .005;


% initializations
pawNames = {'leading fore', 'lagging fore', 'leading hind', 'lagging hind'};
colors = hsv(4);


%% getting distances of penultimate steps for each paw, as well as predicted distances to obs
% also get lengths of penultimate and pre-penultimate steps
[plantingDistances, predictedDistances, penultLength, prePenultLength] = deal(nan(length(data), 4));

for i = 1:length(data)
    
    seq = [data(i).pawOverSequence]'; % sequence with which paws cross obtalce
    if all(ismember(seq(1:2), [2 3])) && all(ismember(seq(3:4), [1 4])) % if first two paws over are forelimb (2,3) and second two paws over are hindlimb (1,4)
        for j = [data(i).pawOverSequence]'
            modSteps = data(i).modStepNum(j);
            
%             if modSteps~=1 % if paw went over in one big step, then penultimate step was complete before obstacle was even reached!
%                 plantingDistances(i,j) = data(i).modifiedLocations{j}(end,1,1); % starting location of step over obstacle
%                 predictedDistances(i,j) = data(i).modifiedLocations{j}(end-1,1,1) + data(i).modPredictedLengths(modSteps-1,j); % take the starting position of the penultimate step, and add the predicted lengt of the penultimate step yo!
%             end

            plantingDistances(i,j) = data(i).modifiedLocations{j}(end,1,1); % starting location of step over obstacle
            if modSteps~=1 % if paw went over in one big step, then penultimate step was complete before obstacle was even reached!    
                predictedDistances(i,j) = data(i).modifiedLocations{j}(end-1,1,1) + data(i).modPredictedLengths(modSteps-1,j); % take the starting position of the penultimate step, and add the predicted lengt of the penultimate step yo!
            else
                predictedDistances(i,j) = data(i).controlLocations{j}(end,1,1) + data(i).controlPredictedLengths(end,j); % take the starting position of the penultimate step, and add the predicted lengt of the penultimate step yo!
            end
            
            % step lengths
            if modSteps==1
                prePenultLength(i,j) = data(i).controlSwingLengths(end-1,j);
                penultLength(i,j) = data(i).controlSwingLengths(end,j);
            elseif modSteps==2
                prePenultLength(i,j) = data(i).controlSwingLengths(end,j);
                penultLength(i,j) = data(i).modifiedSwingLengths(1,j);
            elseif modSteps>2
                prePenultLength(i,j) = data(i).modifiedSwingLengths(end-2,j);
                penultLength(i,j) = data(i).modifiedSwingLengths(end-1,j);
            end
            
        end
    end
end
fprintf('%.2f of trials had a goofy sequence of paws crossing obs...\n', mean(all(isnan(plantingDistances),2)));

%% scatter planting distance vs obs height and velocity

figure('name', figTitle, 'Color', 'white', 'MenuBar', 'none', 'Position', [2000 50 800 300], 'inverthardcopy', 'off')

vars = {[data.obsHeightsVid], [data.vel]};
varNames = {'osbtacle height (mm)', 'velocity (m/s)'};

for var = 1:length(vars)
    
    subplot(1,2,var)
    
    scatters = nan(1,4);
    for i = pawsToPlot
        scatters(i) = scatter(vars{var}, plantingDistances(:,i), 20, colors(i,:), 'filled', 'markerfacealpha', .2); hold on
        validBins = ~isnan(plantingDistances(:,i));
        fit = polyfit(vars{var}(validBins), plantingDistances(validBins,i)', 1);
        plot(vars{var}, polyval(fit, vars{var}), 'linewidth', 4, 'color', colors(i,:));
    end
    
    set(gca, 'YLim', yLims)
    xlabel(varNames{var})
    ylabel('planting distance (m)')
end

legend(scatters(pawsToPlot), pawNames(pawsToPlot))

%% histograms comparing predicted and actual distance to obstacle


figure('name', figTitle, 'Color', 'white', 'MenuBar', 'none', 'Position', [2000 50 500 800], 'inverthardcopy', 'off')

binEdges = xLims(1)-binSize:binSize:xLims(2)+binSize;

for i = 1:4
    subplot(4,1,i)
    
    histogram(plantingDistances(:,i), binEdges, 'FaceColor', colors(i,:)); hold on
    histogram(predictedDistances(:,i), binEdges, 'FaceColor', [.5 .5 .5]); hold on
    
    title(pawNames{i})
    set(gca, 'box', 'off', 'YColor', [1 1 1], 'XLim', xLims + [-.02 .02])
end

xlabel('planting distance to obstacle (m)')



%% histograms comparing penultimate and pre-penultimate steps

figure('name', figTitle, 'Color', 'white', 'MenuBar', 'none', 'Position', [2000 50 500 800], 'inverthardcopy', 'off')

binEdges = 0:binSize:.1;

for i = 1:4
    subplot(4,1,i)
    
    histogram(penultLength(:,i), binEdges, 'FaceColor', colors(i,:)); hold on
    histogram(prePenultLength(:,i), binEdges, 'FaceColor', [.5 .5 .5]); hold on
    
    title(pawNames{i})
    set(gca, 'box', 'off', 'YColor', [1 1 1], 'XLim', [binEdges(1) binEdges(end)])
end

xlabel('step length (m)')





