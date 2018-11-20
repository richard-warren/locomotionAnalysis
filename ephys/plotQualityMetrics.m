function plotQualityMetrics(session)



% SETTINGS
% general
spkWindow = [-.5 1.5]; % ms pre and post spike time to plot
spkNum = 10000; % take this many of all spikes to analyze (to save system resources)
showFigures = 'off';
% highPassFreq = 300;

% snr settings
snrYLims = [0 20];

% probe view settings
xSpacing = 4;
ySpacing = 30;
timeBinNum = 4;
minFiringRate = 0.5; % (hz) time bins with firing rate less than minFr will not display average trace

% spike trace settings
sampleTraceLength = .15;
verticalSpacing = 800;
traceTimeBins = 6;

% firing rate settings
frTimeBinNum = 1000;
binSmoothing = 10;

% false positive settings
fpTimeBinNum = 200;
fpBinSize = 5*60; % s
refractoryPeriod = .003;
censoredPeriod = .001;

% autocorrelogram
xcorrBinWidth = .001;
xcorrWidth = .04;


% COLLECT SESSION DATA
fprintf('%s: collecting data... ', session)
ephysInfo = getSessionEphysInfo(session);
for i = fieldnames(ephysInfo)'; eval([i{1} '=ephysInfo.' i{1} ';']); end % extract field names as variables
meanColor = mean(copper(3),1);
spkWindowInds = int64((spkWindow(1)/1000*fs) : (spkWindow(2)/1000*fs));
load(fullfile(getenv('OBSDATADIR'), 'ephys', 'channelMaps', 'kilosort', [mapFile '.mat']), ...
    'xcoords', 'ycoords')


% get spike times for good units
addpath(fullfile(getenv('GITDIR'), 'npy-matlab'))
addpath(fullfile(getenv('GITDIR'), 'analysis-tools'))
allSpkInds = readNPY(fullfile(getenv('OBSDATADIR'), 'sessions', session, ephysFolder, 'spike_times.npy'));
clusters = readNPY(fullfile(getenv('OBSDATADIR'), 'sessions', session, ephysFolder, 'spike_clusters.npy'));
clusterGroups = readtable(fullfile(getenv('OBSDATADIR'), 'sessions', session, ephysFolder, 'cluster_groups.csv'));
unit_ids = clusterGroups.cluster_id(strcmp(clusterGroups.group, 'good'));


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


for cell = 1:length(unit_ids)

    fprintf('plotting cell %i... ', unit_ids(cell))
    figure('Name', sprintf('%s cell %i', session, unit_ids(cell)), 'Visible', showFigures, ...
        'Color', 'white', 'Position', [2000 45 1601 865]); hold on
    spkInds = allSpkInds(clusters==unit_ids(cell));

    
    % EXTRACT WAVEFORM ACROSS CHANNELS
    spkIndsSub = spkInds(round(linspace(2, length(spkInds), spkNum))); % get subpopulation of spikes evenly spaced out
    spkIndsSubAll = uint64(repmat(spkWindowInds,1,length(spkIndsSub)) + int64(repelem(spkIndsSub, length(spkWindowInds)))');

    allWaveforms = getVoltage(data, 1:channelNum, spkIndsSubAll);
    allWaveforms = reshape(allWaveforms, channelNum, length(spkWindowInds), []);
    allWaveforms = permute(allWaveforms, [3 1 2]);
    allWaveforms = allWaveforms - allWaveforms(:,:,1); % subtract beginning of trace from rest of trace (a hack of a high pass filter, lol)
    
    % find best channel and get voltage for that channel
    meanWaveform = squeeze(mean(allWaveforms,1));
    [~, bestChannel] = max(peak2peak(meanWaveform,2));
    channelData = getVoltage(data, bestChannel, 1:smps);


    % PLOT SPIKE SHAPES ON PROBE
    subplot(4,4,[1 5 9]); hold on
    timeBins = discretize(double(spkIndsSub), timeBinNum);
    colors = copper(timeBinNum);
    sameShankInds = find(abs(xcoords - xcoords(bestChannel))<50);
    
    for j = sameShankInds'
        for i = 1:timeBinNum
            firingRate = sum(timeBins==i) / (range(timeStamps)/timeBinNum);
            if firingRate>minFiringRate % don't plot average trace if rate of spikes in bin is too low, which happens when the unit is lost
                trace = squeeze(mean(allWaveforms(timeBins==i,j,:),1));
                plot(xcoords(j)*xSpacing + spkWindowInds, ...
                    ycoords(j)*ySpacing + trace, ...
                    'Color', colors(i,:), 'LineWidth', 2)
            end
        end
        text(double(xcoords(j)*xSpacing+spkWindowInds(1)), ...
                ycoords(j)*ySpacing, ...
                num2str(j))
    end
    set(gca, 'visible', 'off')



    % PLOT SAMPLE TRACES
    subplot(4,4,14:16)
    colors = copper(traceTimeBins);

    traceTimes = linspace(timeStamps(1), timeStamps(end), traceTimeBins+1);
    for i = 1:traceTimeBins

        % plot raw trace
        traceBins = timeStamps>=traceTimes(i) & timeStamps<traceTimes(i)+sampleTraceLength;
        times = timeStamps(traceBins);
        trace = channelData(traceBins);
        plot(times-times(1), trace + (traceTimeBins-i)*verticalSpacing, 'Color', [.5 .5 .5], 'LineWidth', 1); hold on % plot raw trace

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
    xLims = [0 range(timeStamps)/60];

    amplitudes = peak2peak(squeeze(allWaveforms(:,bestChannel,:)),2);
    stdev = std(channelData);
    minutes = (timeStamps(spkIndsSub)-timeStamps(1)) / 60;
    scatter(minutes, amplitudes/stdev, 20, copper(length(spkIndsSub)));
    set(gca, 'XLim', xLims, 'YLim', snrYLims, 'XColor', get(gcf, 'color'))
    ylabel(sprintf('channel %i SNR', bestChannel))
    xlabel('time (min)')


    % FIRING RATE OVER TIME
    subplot(4,4,6:8)
    [counts, edges] = histcounts(timeStamps(spkInds), frTimeBinNum);
    binWidth = median(diff(edges));
    spkRates = counts / binWidth;
    spkRates = smooth(spkRates, binSmoothing);
    minutes = (edges(1:end-1) + .5*binWidth - timeStamps(1)) / 60; % get bin centers and convert from seconds to minutes
    plot(minutes, spkRates, 'LineWidth', 2, 'Color', meanColor)
    ylabel('firing rate (Hz)')
    set(gca, 'Box', 'off', 'XLim', xLims, 'XColor', get(gcf, 'color'))


    % FALSE POSITIVE RATE
    subplot(4,4,10:12)

    fpRates = nan(1,fpTimeBinNum);
    centers = linspace(0, timeStamps(spkInds(end)), fpTimeBinNum);

    for i = 1:fpTimeBinNum

        binMin = max(0, centers(i)-.5*fpBinSize);
        binMax = min(timeStamps(spkInds(end)), centers(i)+.5*fpBinSize);

        spkBins = timeStamps(spkInds)>=binMin & timeStamps(spkInds)<binMax;
        violations = sum(diff(timeStamps(spkInds(spkBins)))<refractoryPeriod);
        if violations>0
            rhs = 2 * (refractoryPeriod-censoredPeriod) * length(spkInds)^2 / violations / (binMax-binMin); % what the fuck is this math about??? is this really correct???
            fpRates(i) = .5 - .5*sqrt((rhs-4)/rhs);
        else
            fpRates(i) = 0;
        end
    end

    plot((centers-timeStamps(1))/60, fpRates, 'LineWidth', 2, 'Color', meanColor)
    set(gca, 'box', 'off', 'XLim', xLims, 'YLim', [0 .05])
    ylabel('false positive rate')
    xlabel('time (min)')



    % AUTOCORRELOGRAM
    subplot(4,4,13)
    spksBinary = histcounts(timeStamps(spkInds), 0:xcorrBinWidth:timeStamps(end));
    autocorr = xcorr(spksBinary, round(xcorrWidth/2/xcorrBinWidth), 'coeff');
    autocorr(autocorr==1) = nan;
    bar(autocorr, 'BarWidth', 1, 'EdgeColor', 'none', 'FaceColor', mean(copper(3),1))
    set(gca, 'visible', 'off')
    
    
    
    % save figures
    savefig(fullfile(folder, [session 'cell' num2str(unit_ids(cell))]))
    saveas(gcf, fullfile(folder, [session 'cell' num2str(unit_ids(cell)) '.png']))
end
disp('all done!')
warning('on', 'MATLAB:table:ModifiedAndSavedVarnames')


