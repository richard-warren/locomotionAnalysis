function matchedTimes = findMatchedSignalTimes(sig, sigTimes, eventTimes, window, nbest, showPlot)

% given a continuous signal and event times, find other times within the
% signal whether the shape of the signal matches that of the peri-event
% period // use e.g. to find slow-downs that are matched to the slow-downs
% that occur around reward delivery

% finds matches by sliding the average PSTH for the signal across the
% entire signal and taking the nbest moments where the squared error is 
% minimal and there is no overlap with perievent period (defined by window)

% sig:          continuous signal
% sigTimes:     time stamps for continuous signal // MUST BE MONOTONICALLY INCREASING, sorry :(
% eventsTimes:  times of events
% window:       (s) time surrounding events to measure sig, e.g. [-2 2]
% nbest:        take nbest best matched periods

% inits
if ~exist('showPlot', 'var'); showPlot = false; end

% get average response (kernel)
dt = nanmedian(diff(sigTimes));
xk = window(1) : dt : window(2);  % x axis for kernel
X = xk + eventTimes;
responses = interp1(sigTimes, sig, X);
kernel = nanmean(responses, 1);

% mask out peri-event periods
sigMasked = sig;
maskInds = knnsearch(sigTimes', X(:));
sigMasked(maskInds) = nan;

% mean mean squared error between average response kernel and signal for
% all sigTimes
nk = length(kernel);
mse = nan(1, length(sig));
for i = 1:(length(sig)-nk)
    mse(i) = sum((sigMasked(i:i+nk-1) - kernel).^2);
end

% find best matches
[~, matchedTimes] = findpeaks(-mse, sigTimes, 'SortStr', 'descend', 'MinPeakDistance', 5, 'NPeaks', nbest);
matchedTimes = matchedTimes - xk(1);  % shift to the left to compensate for how diffs was computed from left edge of kernel
matchedTimes = matchedTimes(:);


if showPlot
    
    % psth for events and matched events
    figure('color', 'white', 'position', [222.00 655.00 1024.00 510.00]);
    c = lines(1);

    % events
    subplot(2,3,1); hold on
    plot(xk, responses', 'color', [0 0 0 .1])
    plot(xk, mean(responses,1), 'LineWidth', 3, 'color', [0 0 0])
    set(gca, 'xlim', window);
    xlabel('time (s)')
    ylabel('signal')
    yLims = get(gca, 'ylim');
    title('events (kernel)')
    
    % matched events
    subplot(2,3,2); hold on
    matched = interp1(sigTimes, sig, xk + matchedTimes);
    plot(xk, matched', 'color', [c .1])
    plot(xk, mean(matched, 1), 'LineWidth', 3, 'color', c)
    set(gca, 'xlim', window, 'ylim', yLims);
    title('matched events')
    
    % both
    subplot(2,3,3); hold on
    mn = mean(responses, 1);
    plot(xk, mn, 'LineWidth', 3, 'color', [0 0 0])
    stdev = std(responses, 1);
    patch([xk fliplr(xk)], [(-stdev+mn) fliplr(stdev+mn)], [0 0 0], 'FaceAlpha', .1, 'EdgeColor', 'none')
    mn = mean(matched, 1);
    plot(xk, mn, 'LineWidth', 3, 'color', c)
    stdev = std(matched, 1);
    patch([xk fliplr(xk)], [(-stdev+mn) fliplr(stdev+mn)], c, 'FaceAlpha', .1, 'EdgeColor', 'none')
    set(gca, 'ylim', yLims)
    

    % raw signal
    subplot(2,3,4:6); hold on
    plot(sigTimes, sig, 'color', [0 0 0 .1])
    plot(sigTimes, sigMasked, 'color', [0 0 0])
    yLims = get(gca, 'ylim');
    plot(repmat(eventTimes,1,2), yLims, 'color', c)  % plot event times
    set(gca, 'ylim', yLims, 'xlim', [sigTimes(1) sigTimes(end)])

    for i = 1:length(matchedTimes)
        xsub = xk + matchedTimes(i);
        plot(xsub, kernel, 'color', c, 'linewidth', 2)
    end
end

