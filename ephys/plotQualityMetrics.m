function plotQualityMetrics(session)



% SETTINGS
% general
spkWindow = [-.3 .75]; % ms pre and post spike time to plot
spkNum = 10000; % take this many of all spikes to analyze (to save system resources)
showFigures = 'on';
% highPassFreq = 300;

% SNR settings
snrYLims = [0 20];
snrDt = 1; % (s)
snrWindow = 30; % (s)

% probe view settings
xSpacing = 4;
ySpacing = 30;
timeBinNum = 4;
minFiringRate = 0.5; % (hz) time bins with firing rate less than minFr will not display average trace

% spike trace settings
sampleTraceLength = .15;
verticalSpacing = 200;
traceTimeBins = 6;

% firing rate settings
frWindow = 30; % seconds

% false positive settings
fpTimeBinNum = 200;
fpBinSize = 5*60; % s
refractoryPeriod = .002;
% censoredPeriod = .001;

% autocorrelogram
xcorrBinWidth = .001;
xcorrWidth = .04;


% COLLECT SESSION DATA
fprintf('%s: collecting data... ', session)
ephysInfo = getSessionEphysInfo(session);
for i = fieldnames(ephysInfo)'; eval([i{1} '=ephysInfo.' i{1} ';']); end % extract field names as variables
spkWindowInds = int64((spkWindow(1)/1000*fs) : (spkWindow(2)/1000*fs));
load(fullfile(getenv('OBSDATADIR'), 'ephys', 'channelMaps', 'kilosort', [mapFile '.mat']), ...
    'xcoords', 'ycoords', 'channelNum_OpenEphys')

% only for visualizatoon purpose (the waveform plot)
% xcoords(33:64) = 21;
xcoords = xcoords*0.15;
ycoords = ycoords*0.15;


[allSpkInds, unit_ids, bestChannels] = getGoodSpkInds(session); % get spike times for good units
% bestChannels = getBestChannels(session, ephysInfo);

% function to extract voltage from binary file
getVoltage = @(data, channel, inds) ...
    double(data.Data.Data(channel,inds))*bitVolts; % extract voltage from memmapfile, converting to votlage, and only return specific channel

% load data
data = memmapfile(fullfile(getenv('OBSDATADIR'), 'sessions', session, ephysFolder, [fileNameBase '_CHs.dat']), ...
    'Format', {'int16', [channelNum, smps], 'Data'}, 'Writable', false);

% create folder for figures, deleting old one if it already exists
folder = fullfile(getenv('OBSDATADIR'), 'figures', 'ephys', 'qualityMetrics', session);
if isfolder(folder); cmd_rmdir(folder); end % delete folder if already exists
mkdir(folder);


cellColors = hsv(length(unit_ids))*.8;
getColors = @(length, color) interp2(1:3, [1:2]', cat(1,[0 0 0], color), 1:3, linspace(1,2,length)'); % creates gradient from black to color

for c = 1:length(unit_ids)

    fprintf('plotting cell %i... ', unit_ids(c))
    figure('Name', sprintf('%s cell %i', session, unit_ids(c)), 'Visible', showFigures, ...
        'Color', 'white', 'Position', get(0,'ScreenSize')); hold on
    spkInds = allSpkInds{c};

    
    % EXTRACT WAVEFORM ACROSS CHANNELS
    spkIndsSub = spkInds(round(linspace(2, length(spkInds), spkNum))); % get subpopulation of spikes evenly spaced out
    spkIndsSubAll = uint64(repmat(spkWindowInds,1,length(spkIndsSub)) + int64(repelem(spkIndsSub, length(spkWindowInds)))');

    allWaveforms = getVoltage(data, 1:channelNum, spkIndsSubAll);
    allWaveforms = reshape(allWaveforms, channelNum, length(spkWindowInds), []);
    allWaveforms = permute(allWaveforms, [3 1 2]);
%     allWaveforms = allWaveforms - allWaveforms(:,:,1); % subtract beginning of trace from rest of trace (a hack of a high pass filter, lol)
    
    % find best channel and get voltage for that channel
%     meanWaveform = squeeze(mean(allWaveforms,1));
%     [~, bestChannel] = max(peak2peak(meanWaveform,2));
    channelData = getVoltage(data, bestChannels(c), 1:smps);
%     templates = readNPY(fullfile(getenv('OBSDATADIR'), 'sessions', session, ephysFolder, 'templates.npy'));
%     spike_templates = readNPY(fullfile(getenv('OBSDATADIR'), 'sessions', session, ephysFolder, 'spike_templates.npy'));
    


    % PLOT SPIKE SHAPES ON PROBE
    subplot(4,4,[1 5 9]); hold on
    timeBins = discretize(double(spkIndsSub), timeBinNum);
    colors = getColors(timeBinNum, cellColors(c,:));
    sameShankInds = find(abs(xcoords - xcoords(bestChannels(c)))<50);
    
    for j = 1:64
        for i = 1:timeBinNum
            firingRate = sum(timeBins==i) / (range(timeStamps)/timeBinNum);
            if firingRate>minFiringRate % don't plot average trace if rate of spikes in bin is too low, which happens when the unit is lost
                trace = squeeze(mean(allWaveforms(timeBins==i,j,:),1));
                plot(xcoords(j)*xSpacing + double(spkWindowInds), ...
                    ycoords(j)*ySpacing + trace*0.3, ...
                    'Color', colors(i,:), 'LineWidth', 2)
            end
        end
        if j==bestChannels(c); textColor='red'; else; textColor='black'; end
        text(double(xcoords(j)*xSpacing+spkWindowInds(1)), ...
                ycoords(j)*ySpacing, ...
                num2str(find(channelNum_OpenEphys == j)), 'Color', textColor)
    end
    set(gca, 'visible', 'off')


    % PLOT SAMPLE TRACES
    subplot(4,4,14:16)
    colors = getColors(traceTimeBins, cellColors(c,:));

    traceTimes = linspace(timeStamps(1), timeStamps(end), traceTimeBins+1);
    for i = 1:traceTimeBins

        % plot raw trace
        traceBins = timeStamps>=traceTimes(i) & timeStamps<traceTimes(i)+sampleTraceLength;
        times = timeStamps(traceBins);
        trace = channelData(traceBins);
        plot(times-times(1), trace + (traceTimeBins-i)*verticalSpacing, 'Color', [.5 .5 .5], 'LineWidth', 0.5); hold on % plot raw trace

        % overlay spikes
        traceSpkInds = spkInds(timeStamps(spkInds)>=traceTimes(i) & timeStamps(spkInds)<traceTimes(i)+sampleTraceLength);
        spkBins = repmat(spkWindowInds,length(traceSpkInds),1) + int64(traceSpkInds); % matrix where each row is inds for given spike
        for j = 1:length(traceSpkInds)
            plot(timeStamps(spkBins(j,:))-times(1), channelData(spkBins(j,:)) + (traceTimeBins-i)*verticalSpacing, ...
                'Color', colors(i,:), 'LineWidth', 2); % overlay spikes! lol
        end
    end
    set(gca, 'visible', 'off', 'XLim', [0 range(times)])




    % AMP OVER TIME
    subplot(4,4,2:4)
    xLims = [min(timeStamps) max(timeStamps)]/60;

    amplitudes = peak2peak(squeeze(allWaveforms(:,bestChannels(c),:)),2);
    stdev = std(channelData);
    
    scatter(timeStamps(spkIndsSub)/60, amplitudes/stdev, ...
        20, getColors(length(spkIndsSub), cellColors(c,:)), 'MarkerEdgeAlpha', .4);
    set(gca, 'XLim', xLims, 'YLim', snrYLims, 'XColor', get(gcf, 'color'))
    ylabel(sprintf('channel %i SNR', bestChannels(c)))
    xlabel('time (min)')


    % FIRING RATE OVER TIME
    subplot(4,4,6:8)
    [counts, edges] = histcounts(timeStamps(spkInds), min(timeStamps):1:max(timeStamps));
    binWidth = median(diff(edges));
    spkRates = counts / binWidth;
    spkRates = smooth(spkRates, frWindow/binWidth);
    minutes = (edges(1:end-1) + .5*binWidth) / 60; % get bin centers and convert from seconds to minutes
    plot(minutes, spkRates, 'LineWidth', 2, 'Color', cellColors(c,:))
    ylabel('firing rate (Hz)')
    set(gca, 'Box', 'off', 'XLim', xLims, 'XColor', get(gcf, 'color'))


    % FALSE POSITIVE RATE
    subplot(4,4,10:12)
    yyaxis right
    fpRates = nan(1,fpTimeBinNum);
    centers = linspace(0, timeStamps(spkInds(end)), fpTimeBinNum);

    for i = 1:fpTimeBinNum

        binMin = max(0, centers(i)-.5*fpBinSize);
        binMax = min(timeStamps(spkInds(end)), centers(i)+.5*fpBinSize);

        spkBins = timeStamps(spkInds)>=binMin & timeStamps(spkInds)<binMax;
        violations = sum(diff(timeStamps(spkInds(spkBins)))<refractoryPeriod);
%         if violations>0
%             rhs = 2 * (refractoryPeriod-censoredPeriod) * sum(spkBins)^2 / violations / (binMax-binMin); % what the fuck is this math about??? is this really correct???
%             fpRates(i) = .5 - .5*sqrt((rhs-4)/rhs);
%         else
%             fpRates(i) = 0;
%         end
        fpRates(i) = violations/sum(spkBins);
    end

    plot(centers/60, fpRates, 'LineWidth', 2, 'Color', 'black')
    set(gca, 'box', 'off', 'XLim', xLims, 'YLim', [0 .05], 'YColor', 'black');
    ylabel('false positive rate')
    xlabel('time (min)')
    
    
    % 5th percentile SNR
    times = min(timeStamps):snrDt:max(timeStamps);
    percentiles = nan(1,length(times));
    for i = 1:length(times)
        spkBins = timeStamps(spkIndsSub)>times(i)-snrWindow*.5 & timeStamps(spkIndsSub)<=times(i)+snrWindow*.5;
        binSnrs = amplitudes(spkBins) / stdev; % !!! change to bin std, not global std?
        if sum(spkBins)/snrWindow > minFiringRate
            percentiles(i) = prctile(binSnrs, 5);
        end
    end
    yyaxis left
    plot(times/60, percentiles, 'LineWidth', 2, 'Color', cellColors(c,:))
    ylabel('5th percentile SNR')
    set(gca, 'YColor', cellColors(c,:), 'YLim', snrYLims*.5, 'Box', 'off');



    % AUTOCORRELOGRAM
    subplot(4,4,13)
    spksBinary = histcounts(timeStamps(spkInds), 0:xcorrBinWidth:timeStamps(end));
    autocorr = xcorr(spksBinary, round(xcorrWidth/2/xcorrBinWidth), 'coeff');
    autocorr(autocorr==1) = nan;
    bar(autocorr, 'BarWidth', 1, 'EdgeColor', 'none', 'FaceColor', cellColors(c,:))
    set(gca, 'visible', 'off')
    
    
    % save figures
    blackenFig
    savefig(fullfile(folder, [session 'cell' num2str(unit_ids(c))]))
    saveas(gcf, fullfile(folder, [session 'cell' num2str(unit_ids(c)) '.png']))
end

% plotting all channels
showChannelsOverTime(session, 5, true, fullfile(folder, [session 'allUnits.png']), bestChannels);


% check if these unit_ids have already been written to cellData.csv file
fileName = fullfile(getenv('OBSDATADIR'), 'sessions', session, 'cellData.csv');
alreadyWritten = false;
if exist(fileName, 'file')
    csvTable = readtable(fileName);
    alreadyWritten = isequal(csvTable.unit_id, unit_ids);
end

% make new cellData.csv file if it doesn't already exist
if ~alreadyWritten
    disp('writing cellData.csv metadata file...')
    csvTable = cell2table(cell(0,6), ...
        'VariableNames', {'unit_id', 'include', 'timeStart', 'timeEnd', 'location', 'notes'});
    warning('off', 'MATLAB:table:RowsAddedExistingVars')
    csvTable(1:length(unit_ids),1) = table(unit_ids);
    writetable(csvTable, fileName)
    warning('on', 'MATLAB:table:RowsAddedExistingVars')
end

disp('all done!')


