function [plotInd_new] = plotContinuousData(data, unitSpkRate, plotInd, bins, frameTimeStamps, neuralTimes, opts)

% settings
s.rows = 4;
s.cols = 4;
s.xLabel = ' ';
s.yLabel = ' ';
s.dataName = ' ';
s.xcorrFrequency = 250; % Hz
s.scatterSampleNumber = 1000;
s.scatSize = 20; % scatter size
s.scatAlpha = .5; % scatter circle transparency
s.lineAlpha = 1; % line transparency
s.lineColor = [1, .4, .4]; % scatter point color
s.pointColor = [.08, .6, 1]; % fit line color
s.rawTraceNumber = 5; % number of raw traces plotted
s.rawTraceDuration = 4; % second
s.fillMissingMethod = 'spline'; % method used for filling any nans in the data
s.zscoreThresh = 10; % zscore thresh to get rid of outliers

% reassign settings contained in opts
if exist('opts', 'var'); for i = 1:2:length(opts); s.(opts{i}) = opts{i+1}; end; end



if exist('bins', 'var')
    
    binnedData = data(bins);
    binnedDataZscores = zscore(binnedData);
    
    % plot distributions of the data and zscores of the data
    % disp('plotting sub fig 1...');
    subplot(s.rows, s.cols, plotInd);
    
    currentAxes(1) = subplot(s.rows, s.cols, plotInd);
    tempPosition = [currentAxes(1).Position];
    histogram(binnedData(abs(binnedDataZscores) < 10));
    xlabel(s.dataName);
    ylabel('counts');
    title(s.dataName);
    box off
    
    currentAxes(2) = axes('Position', [tempPosition(1) + 0.115, tempPosition(2) + 0.11, tempPosition(3)*.2, tempPosition(4)*.2]);
    box on
    plot(binnedDataZscores(abs(binnedDataZscores) < 10));
    ylabel('zscore');
    xlabel('inds');
    axis tight
    
    clear currentAxes
    
    
    % plot cross correlations between the data and neural firing rates
    % disp('plotting sub fig 2...');
    plotInd = plotInd + 1;
    subplot(s.rows, s.cols, plotInd);
    binnedData = data(bins);
    
    % cross correlation sanity check
    if any(isnan(data))
        warning('Continuous data contains NAN!!!')
        [data,data_filled] = fillmissing(data, s.fillMissingMethod);
        fprintf('Dealt with fillmissing function. Method = %s', s.fillMissingMethod);
    end
    
    if any(isnan(unitSpkRate))
        warning('Neural data contains NAN!!!')
        [unitSpkRate, unitSpkRate_filled] = fillmissing(data, s.fillMissingMethod);
        fprintf('Dealt with fillmissing function. Method = %s', s.fillMissingMethod);
    end    
    
    % start data processing for cross correlation 
    startInds = find(diff(bins) == 1) + 1;
    endInds = find(diff(bins) == -1);
    
    while(endInds(1) < startInds(1))
        endInds(1) = [];
    end
    
    while(startInds(end) > endInds(end))
        startInds(end) = [];
    end

    % calculate chunk data and corresponding neural firing rates
    if length(startInds) ~= length(endInds)
        warning('Cross correlation: chunk data startinds and endInds DO NOT match!!')
    else
        FRInterp = [];
        for i = 1:length(startInds)
            k = 1;
            chunkData = data(startInds(i):endInds(i));          
            startTime = frameTimeStamps(startInds(i));
            endTime = frameTimeStamps(endInds(i));
            chunkFR = unitSpkRate(neuralTimes > startTime & neuralTimes < endTime);
            chunkFRInterp = interp1(1:length(chunkFR), chunkFR, linspace(1, length(chunkFR), endInds(i) - startInds(i) + 1));
            FRInterp = [FRInterp, chunkFRInterp];
            
            if any(abs(zscore(chunkData)) >= 10)
                continue
            end           
            
            if length(chunkFRInterp) ~= length(chunkData)
                warning('Cross correlation: length of chunk FR DOES NOT match length of data!!!')
            else
                [c(k, :), lags] = xcorr(chunkData-mean(chunkData), chunkFRInterp-mean(chunkFRInterp), 2*s.xcorrFrequency, 'normalized');
                k = k+1;
            end
        end
        
        % plot!!!
        x = lags/s.xcorrFrequency;
        plot(x, nanmean(c, 1), '-', 'LineWidth', 4);
        xlabel('time (second)');
        ylabel('mean corr coef');
        box off
        axis tight;
      
        yLimit=get(gca,'ylim');
        xLimit=get(gca,'xlim');
        yRange = yLimit(2) - yLimit(1);
        
        maxind = find(nanmean(c, 1) == max(nanmean(c, 1)));        
        line([x(maxind), x(maxind)], yLimit, 'linewidth', 2, 'color', s.lineColor);
        minind = find(nanmean(c, 1) == min(nanmean(c, 1)));      
        line([x(minind), x(minind)], yLimit, 'linewidth', 2, 'color', [1 0.77 0.16]);
        
        text(xLimit(1), yLimit(2), ['maximum.x = ', num2str(x(maxind))], 'Color', s.lineColor, 'FontWeight', 'bold');
        text(xLimit(1), yLimit(2)-yRange*0.1, ['corrcoef = ', num2str(max(mean(c, 1)))], 'Color', s.lineColor, 'FontWeight', 'bold');
        text(xLimit(1), yLimit(2)-yRange*0.2, ['minimum.x = ', num2str(x(minind))], 'Color', [1 0.77 0.16], 'FontWeight', 'bold');
        text(xLimit(1), yLimit(2)-yRange*0.3, ['corrcoef = ', num2str(min(mean(c, 1)))], 'Color', [1 0.77 0.16], 'FontWeight', 'bold');    
    end
    
    
    
    % scatters b/w the data and neural firing rates
    % disp('plotting sub fig 3...');
    plotInd = plotInd + 1;
    subplot(s.rows, s.cols, plotInd);
    binnedData = data(bins);
    inds = find(abs(zscore(binnedData)) >= 10);
    binnedData(inds) = [];
    FRInterp(inds) = [];
   
    if length(binnedData) ~= length(FRInterp)
        warning('Length of binned data DOES NOT match length of interpolated neural firing rates.');
    else               
        samples = datasample(1:length(binnedData), s.scatterSampleNumber);
        
        scatter(binnedData(samples), FRInterp(samples), s.scatSize, s.pointColor, 'filled', ...
            'MarkerEdgeAlpha', s.scatAlpha, 'MarkerFaceAlpha', s.scatAlpha);
        hold on; box off;
        xlabel(s.dataName);
        ylabel('neural firing rate (spk/s)')
        ylim([0, max(FRInterp(samples))]);
        
        if size(binnedData, 2) ~= 1; binnedData = binnedData'; end
        if size(FRInterp, 2) ~= 1; FRInterp = FRInterp'; end
        fit = polyfit(binnedData, FRInterp, 1);
        plot(binnedData, polyval(fit, binnedData), 'linewidth', 4, 'color', [s.lineColor s.lineAlpha]);
        corrs = corr(binnedData, FRInterp);
        slopes = fit(1);
        
        yLimit=get(gca,'ylim');
        xLimit=get(gca,'xlim');
        yRange = yLimit(2) - yLimit(1);
        
        text(xLimit(1), yLimit(2), ['R^2 = ', num2str(corrs^2)]);
        text(xLimit(1), yLimit(2)-yRange*0.1, ['slope = ', num2str(slopes)]);
    end
    
    
    % plotting raw traces for data with the most zscores
    % disp('plotting sub fig 4...');
    plotInd = plotInd + 1;
    subplot(s.rows, s.cols, plotInd);hold on
    box off
    
    displacement = 0;
    samples = [];
    samples = datasample(1:length(startInds), s.rawTraceNumber);
    
    for i = 1:s.rawTraceNumber
        
        % get interpolated raw traces for data and neural firing rates
        chunkData = data(startInds(samples(i)):endInds(samples(i)));
        chunkData = interp1(1:length(chunkData), chunkData, linspace(1, length(chunkData), 200));
        
        startTime = frameTimeStamps(startInds(samples(i)));
        endTime = frameTimeStamps(endInds(samples(i)));
        chunkFR = unitSpkRate(neuralTimes > startTime & neuralTimes < endTime);
        chunkFRInterp = interp1(1:length(chunkFR), chunkFR, linspace(1, length(chunkFR), 200));
        
        % plot!!!
        chunkDataRawTrace = chunkData/(max(chunkData) - min(chunkData)) - mean(chunkData/(max(chunkData) - min(chunkData)));
        plot(chunkDataRawTrace + displacement, '-', 'linewidth', 1, 'color', [s.pointColor]);
        hold on
        chunkFRInterpRawTrace = chunkFRInterp/(max(chunkFRInterp) - min(chunkFRInterp)) - mean(chunkFRInterp/(max(chunkFRInterp) - min(chunkFRInterp)));
        plot(chunkFRInterpRawTrace + displacement, '-', 'linewidth', 1, 'color', [s.lineColor]);
        axis tight
        
        displacement = displacement - 1.5;
    end
    
    xlabel('interpolated inds');
    legend(s.dataName, 'neural firing rates');
    
    plotInd_new = plotInd + 1;
  
    
    
else  
    % plot distributions of the data and zscores of the data
    disp('plotting sub fig 1...');
    subplot(s.rows, s.cols, plotInd);
    
    currentAxes(1) = subplot(s.rows, s.cols, plotInd);
    tempPosition = [currentAxes(1).Position];
    histogram(data(bins));
    xlabel(s.dataName);
    ylabel('counts');
    title(s.dataName);
    box off
    
    currentAxes(2) = axes('Position', [tempPosition(1) + 0.115, tempPosition(2) + 0.11, tempPosition(3)*.2, tempPosition(4)*.2]);
    box on
    x = (1:length(data))/s.xcorrFrequency;
    plot(x, zscore(data));
    ylabel('zscore');
    xlabel('frame');
    axis tight
    
    clear currentAxes
    
    
    % plot cross correlations between the data and neural firing rates
    disp('plotting sub fig 2...');
    plotInd = plotInd + 1;
    currentAxes = subplot(s.rows, s.cols, plotInd);
    if length(data) == length(unitSpkRate)
        
        % sanity check
        if any(isnan(data))
            warning('Continuous data contains NAN!!!')
            [data,data_filled] = fillmissing(data, s.fillMissingMethod);
            fprintf('Dealt with fillmissing function. Method = %s', s.fillMissingMethod);
        end
        
        if any(isnan(unitSpkRate))
            warning('Neural data contains NAN!!!')
            [unitSpkRate, unitSpkRate_filled] = fillmissing(data, s.fillMissingMethod);
            fprintf('Dealt with fillmissing function. Method = %s', s.fillMissingMethod);
        end
        
        % plot!
        if size(data, 2) ~= 1; data = data'; end
        if size(unitSpkRate, 2) ~= 1; unitSpkRate = unitSpkRate'; end
        [c, lags] = xcorr(data-mean(data), unitSpkRate-mean(unitSpkRate), 2*s.xcorrFrequency, 'normalized');
        plot(lags/s.xcorrFrequency, c, '-', 'LineWidth', 4);
        xlabel('time (second)');
        ylabel('corr coef');
        box off
        axis tight;
        
        
        yLimit=get(gca,'ylim');
        xLimit=get(gca,'xlim');
        yRange = yLimit(2) - yLimit(1);
        
        maxind = find(c == max(c));
        x = lags/s.xcorrFrequency;
        line([x(maxind), x(maxind)], yLimit, 'linewidth', 2, 'color', s.lineColor);
        minind = find(c == min(c));
        line([x(minind), x(minind)], yLimit, 'linewidth', 2, 'color', [1 0.77 0.16]);
        
        text(xLimit(1), yLimit(2), ['maximum.x = ', num2str(x(maxind))], 'Color', s.lineColor, 'FontWeight', 'bold');
        text(xLimit(1), yLimit(2)-yRange*0.1, ['corrcoef = ', num2str(max(c))], 'Color', s.lineColor, 'FontWeight', 'bold');
        text(xLimit(1), yLimit(2)-yRange*0.2, ['minimum.x = ', num2str(x(minind))], 'Color', [1 0.77 0.16], 'FontWeight', 'bold');
        text(xLimit(1), yLimit(2)-yRange*0.3, ['corrcoef = ', num2str(min(c))], 'Color', [1 0.77 0.16], 'FontWeight', 'bold');
        
    else
        display('length of the data DOES NOT match length of unit firing rate!');
    end
    
    
    
    % scatters b/w the data and neural firing rates
    disp('plotting sub fig 3...');
    plotInd = plotInd + 1;
    currentAxes = subplot(s.rows, s.cols, plotInd);
    inds = datasample(1:length(data), s.scatterSampleNumber);
    
    scatters = scatter(data(inds), unitSpkRate(inds), s.scatSize, s.pointColor, 'filled', ...
        'MarkerEdgeAlpha', s.scatAlpha, 'MarkerFaceAlpha', s.scatAlpha);
    hold on; box off;
    xlabel(s.dataName);
    ylabel('neural firing rate (spk/s)')
    ylim([0, max(unitSpkRate(inds))]);
    
    fit = polyfit(data, unitSpkRate, 1);
    lines = plot(data, polyval(fit, data), 'linewidth', 4, 'color', [s.lineColor s.lineAlpha]);
    corrs = corr(data, unitSpkRate);
    slopes = fit(1);
    
    yLimit=get(gca,'ylim');
    xLimit=get(gca,'xlim');
    yRange = yLimit(2) - yLimit(1);
    
    text(xLimit(1), yLimit(2), ['R^2 = ', num2str(corrs^2)]);
    text(xLimit(1), yLimit(2)-yRange*0.1, ['slope = ', num2str(slopes)]);
    
    
    % plotting raw traces for data with the most zscores
    disp('plotting sub fig 4...');
    plotInd = plotInd + 1;
    subplot(s.rows, s.cols, plotInd);hold on
    box off
    
    [~, inds] = sort(zscore(data), 'descend');
    lineAlpha = linspace(1, 0.1, s.rawTraceNumber);
    indWindow = s.rawTraceDuration*s.xcorrFrequency/2;
    colors = hsv(s.rawTraceNumber);
    displacement = 0;
    
    for i = 1:s.rawTraceNumber
        
        x = linspace(-s.rawTraceDuration/2, s.rawTraceDuration/2, s.rawTraceDuration*s.xcorrFrequency+1);
        if (inds(1) < s.rawTraceDuration*s.xcorrFrequency/2) || (inds(i) > (length(data) - s.rawTraceDuration*s.xcorrFrequency/2))
            continue
        else
            plot(x, data(inds(1)-indWindow : inds(1)+indWindow) + displacement, '-', 'linewidth', 2, 'color', [s.pointColor, lineAlpha(i)]);hold on
        end
        
        inds(inds >= inds(1)-indWindow & inds <= inds(1)+indWindow) = [];
        displacement = displacement - max(data(inds(1)-indWindow : inds(1)+indWindow)) + min(data(inds(1)-indWindow : inds(1)+indWindow));
    end
    
    xlabel('time (second)');
    ylabel(s.dataName);
    
    plotInd_new = plotInd + 1;
end

end