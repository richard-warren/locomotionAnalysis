function plotSensoryDependenceSuccess(sessionInfo, data)
% plot sensory dependence of obstacle avoidance

% settings
touchThresh = 1; % if paw contacts obs for more than touchThresh frames, trial is considered touching

% initializations
sessionInfo = readtable([getenv('OBSDATADIR') 'sessions\sessionInfo.xlsx'], 'Sheet', 'whiskerTrimNotes');
sessionInfo = sessionInfo(sessionInfo.include==1 & ~cellfun(@isempty, sessionInfo.session),:);





% PLOT!

% settings
lineEdges = [-.25 .25];
lineWidth = 3;
colors = winter(2);
circSize = 80;
successYLims = [0 1];
velYLim = [.2 .6];
jitterRange = .25;

% determine whether trial was successful
isSuccess = num2cell((([data.frontTouchingFrames] + [data.topTouchingFrames]) < touchThresh) & ~[data.isWheelBreak]);
[data.isSuccess] = isSuccess{:};


figure('color', 'white', 'menubar', 'none', 'inverthardcopy', 'off', 'position', [100 100 1280/2 720]);


mice = unique({data.mouse});

for wisk = [false true]
    for light = [false true]
        
        successRates = nan(1,length(mice));
        vels = nan(1,length(mice));
        jitters = linspace(-jitterRange/2, jitterRange/2, length(mice));
        
        for mouse = 1:length(mice)
            
            % get session performance
            trialBins = [data.hasWiskers]==wisk & ...
                ([data.isLightOn]==light) & ...
                strcmp(mice{mouse}, {data.mouse}) & ...
                [data.trial]<=maxTrialNum;
            trialBinsNoBreak = trialBins & ~[data.isWheelBreak];
            successRates(mouse) = sum([data(trialBins).isSuccess]) / sum(trialBins);
            vels(mouse) = mean([data(trialBinsNoBreak).vel]);
            
            % plot
            subplot(2,2,wisk+1)
            scatter(light-.5 + jitters(mouse), successRates(mouse), circSize, colors(light+1,:), 'filled'); hold on
            subplot(2,2,wisk+3)
            scatter(light-.5 + jitters(mouse), vels(mouse), circSize, colors(light+1,:), 'filled'); hold on
        end
        
        % plot means
        subplot(2,2,wisk+1)
        line(lineEdges+light-.5, [mean(successRates) mean(successRates)], 'linewidth', lineWidth, 'color', 'black')
        subplot(2,2,wisk+3)
        line(lineEdges+light-.5, [mean(vels) mean(vels)], 'lineWidth', lineWidth, 'color', 'black')
    end
end

% pimp figs
subplot(2,2,1)
set(gca, 'xcolor', 'white', 'ylim', successYLims)
ylabel('success rate')

subplot(2,2,2)
set(gca, 'xcolor', 'white', 'ylim', successYLims)

subplot(2,2,3)
set(gca, 'xtick', [-.5 .5], 'xticklabel', {'no light', 'light'}, 'ylim', velYLim)
xlabel('no whiskers')
ylabel('velocity (m/s)')

subplot(2,2,4)
set(gca, 'xtick', [-.5 .5], 'xticklabel', {'no light', 'light'}, 'ylim', velYLim)
xlabel('whiskers')

saveas(gcf, [getenv('OBSDATADIR') 'figures\sensoryDependennce.png']);















