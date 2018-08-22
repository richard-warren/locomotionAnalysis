function sensoryDependenceSuccessAndSpeed(data)
% plot sensory dependence of obstacle avoidance

% settings
touchThresh = 1; % if paw contacts obs for more than touchThresh frames, trial is considered touching
circSize = 75;
lineThickness = 4;
jitterRange = .2;
lineEdges = [-.2 .2];
dvs = {'success rate', 'speed (m/s)'};
yLims = [0 1; .1 .8];
colorFading = 0.6;



% initializations
totalTouches = cellfun(@mean, {data.totalTouchFramesPerPaw});
isSuccess = totalTouches<touchThresh;

conditionNames = {'whiskers + vision', 'whiskers', 'vision', 'neither'};
wiskBins = strcmp({data.preOrPost}, 'pre');
lightOnBins = [data.isLightOn];
dvBins = {wiskBins & lightOnBins
          wiskBins & ~lightOnBins
          ~wiskBins & lightOnBins
          ~wiskBins & ~lightOnBins};

close all; figure('color', 'white', 'menubar', 'none', 'inverthardcopy', 'off', 'position', [100 100 500 720]);


mice = unique({data.mouse});
colors = hsv(length(conditionNames));
jitters = linspace(-jitterRange/2, jitterRange/2, length(mice));
successes = nan(length(mice), length(conditionNames));
vels = nan(length(mice), length(conditionNames));



for i = 1:length(conditionNames)
    
    clr1 = colors(i,:)+colorFading; clr1(clr1>1)=1;
    clr2 = colors(i,:)-.5*colorFading; clr2(clr2<0)=0;
    mouseColors = interp2(1:3, 1:2, cat(1,clr1,clr2), 1:3, linspace(1,2,length(mice))');
    
    for mouse = 1:length(mice)

        % get session performance
        mouseBins = strcmp({data.mouse}, mice{mouse});
        
        successBins = mouseBins & dvBins{i};
        successes(mouse,i) = mean(isSuccess(successBins));
        
        velBins = ~[data.isWheelBreak] & dvBins{i} & mouseBins;
        vels(mouse,i) = mean([data(velBins).avgVel]);

        % plot mouse avoidance and speed
        subplot(2,1,1)
        scatter(i+jitters(mouse), successes(mouse,i), circSize, mouseColors(mouse,:), 'filled'); hold on
        subplot(2,1,2)
        scatter(i+jitters(mouse), vels(mouse,i), circSize, mouseColors(mouse,:), 'filled'); hold on
    end

    % plot mean avoidance and speed
    subplot(2,1,1)
    line(lineEdges+i, repmat(nanmean(successes(:,i)),1,2), 'linewidth', lineThickness, 'color', colors(i,:))
    
    subplot(2,1,2)
    line(lineEdges+i, repmat(nanmean(vels(:,i)),1,2), 'lineWidth', lineThickness, 'color', colors(i,:))
end

% add mouse lines
for mouse = 1:length(mice)
    
    % success
    subplot(2,1,1)
    ln = line([1:length(conditionNames)] + jitters(mouse), successes(mouse,:), 'color', [.2 .2 .2]);
    uistack(ln, 'bottom')
    % vel
    subplot(2,1,2)
    ln = line([1:length(conditionNames)] + jitters(mouse), vels(mouse,:), 'color', [.2 .2 .2]);
    uistack(ln, 'bottom')
end


% pimp figs
for i = 1:2
    subplot(2,1,i)
    ylabel(dvs{i})
    set(gca, 'XTick', 1:length(conditionNames), 'XTickLabel', conditionNames, 'YLim', yLims(i,:))
end

blackenFig












