function plotPredictorResponse(session, unit, predictor, varargin)

% plot a single predictor for a single session and unit // this is a
% wrapper for plotMovingAvg (for continuous variables) and plotPSTH (for
% event and epoch variables)

% settings
% --------

% x axis
s.eventLims = [-.25 .5];  % (s) x limits for events
s.epochLims = [-.25 1.25];  % (fraction of epoch) x limits for epochs
s.traceLims = [-10 2];  % (s) time pre and post reward to show raw traces for continuous variables  

% instantaneous firing rate kernel
s.kernelRise = .005;     % (s) rise for double exponential kernel
s.kernelFall = .02;      % (s) fall for double exponential kernel
s.kernelSig = .02;       % (s) if a gaussian kernel is used
s.kernel = 'doubleExp';  % 'gauss', or 'doubleExp'

% events (show for continuous vars only)
s.rewardColor  = lines(1);
s.whiskerColor = [174 94 203]/255;

% other
s.maxEpochs = 2000;  % if more than maxEpochs, limit to central maxEpochs
s.color = lines(1);
s.predictorColor = [.2 .2 .2];



% initializations
% ---------------
if exist('varargin', 'var'); for i = 1:2:length(varargin); s.(varargin{i}) = varargin{i+1}; end; end % reassign settings passed in varargin
load(fullfile(getenv('SSD'), 'modelling', 'predictors', [session '_predictors.mat']), 'predictors');
load(fullfile(getenv('OBSDATADIR'), 'sessions', session, 'neuralData.mat'), 'spkTimes', 'unit_ids');
spkTimes = spkTimes{unit_ids==unit};
args = {'xlabel', predictor, 'ylabel', 'firing rate', 'color', s.color};


% plot!
if predictors{predictor, 'type'}=='continuous'
    
    figure('color', 'white', 'position', [177.00 410.00 1563.00 282.00]);
    
    % moving average
    % --------------
    subplot(1,3,1)
    [spkRate, spkRatesTimes] = getFiringRate(spkTimes, 'fs', 2000, ...
        'kernel', s.kernel, 'kernelRise', s.kernelRise, 'kernelFall', s.kernelFall, 'sig', s.kernelSig);
    t = predictors{predictor, 't'}{1};
    spkRate = interp1(spkRatesTimes, spkRate, t);  % put on same time grid
    predictorData = predictors{predictor, 'data'}{1};
    
    yLims = [0 nanmean(spkRate)+2*nanstd(spkRate)];
    plotMovingAvg(predictorData, spkRate, args{:}, 'yLims', yLims, ...
        'plotDensity', false, 'newFig', false, 'color', s.color, 'showScatter', false);
    
    xy = rmmissing([predictorData', spkRate']);
    r = corr(xy(:,1), xy(:,2));
    text(1, 1, sprintf('r = %.2f', r), 'Units', 'normalized', ...
        'HorizontalAlignment', 'right', 'VerticalAlignment', 'top')
    
    % sample traces
    % -------------
    subplot(1,3,2:3); hold on
    rewardTimes = predictors{'reward_all', 'data'}{1};
    whiskerTimes = predictors{'whiskerContact', 'data'}{1};
    
    % determine x axis
    tLims = [t(find(~isnan(spkRate),1,'first')) t(find(~isnan(spkRate),1,'last'))];  % valid times for the unit
    rewardTimes = rewardTimes(rewardTimes>tLims(1) & rewardTimes<tLims(2));
    rewardInd = randsample(5:length(rewardTimes)-5, 1);
    xLims = rewardTimes(rewardInd) + s.traceLims;  % two seconds before and after a pair of rewards
    bins = t>=xLims(1) & t<=xLims(2);
    
    yyaxis left; plot(t(bins), predictorData(bins), 'LineWidth', 1.5, 'color', [s.predictorColor .8])
    ylabel(predictor)
    set(gca, 'YColor', s.predictorColor, 'TickDir', 'out')
    
    yyaxis right; plot(t(bins), spkRate(bins), 'LineWidth', 1.5, 'color', [s.color .8])
    ylabel('firing rate')
    set(gca, 'YColor', s.color, 'TickDir', 'out')
    
    % add lines for events
    yLims = get(gca, 'YLim');
    plot([whiskerTimes whiskerTimes], yLims, '-', 'color', s.whiskerColor)
    scatter(whiskerTimes, repelem(yLims(2), length(whiskerTimes)), 20, s.whiskerColor, 'filled')
    plot([rewardTimes rewardTimes], yLims, '-', 'color', s.rewardColor)
    scatter(rewardTimes, repelem(yLims(2), length(rewardTimes)), 20, s.rewardColor, 'filled')
    
    set(gca, 'xlim', xLims, 'ylim', yLims, 'box', 'off')
    xlabel('time (s)')
    
else
    % PSTH with raster
    % ----------------
    events = predictors{predictor, 'data'}{1};
    plotPSTH(spkTimes, events, 'removeNoSpikeTrials', true, args{:}, 'maxEpochs', s.maxEpochs, ...
        'eventLims', s.eventLims, 'epochLims', s.epochLims, ...
        'kernel', s.kernel, 'kernelRise', s.kernelRise, 'kernelFall', s.kernelFall, 'kernelSig', s.kernelSig);
end



