c = 1

fprintf('plotting cell %i... ', unit_ids(c))
figure('Name', sprintf('%s cell %i', session, unit_ids(c)), 'Visible', showFigures, ...
    'Color', 'white'); hold on
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

for j = sameShankInds
    for i = 1:timeBinNum
        firingRate = sum(timeBins==i) / (range(timeStamps)/timeBinNum);
        if firingRate>minFiringRate % don't plot average trace if rate of spikes in bin is too low, which happens when the unit is lost
            trace = squeeze(mean(allWaveforms(timeBins==i,j,:),1));
            plot(xcoords(j)*xSpacing + double(spkWindowInds), ...
                ycoords(j)*ySpacing + trace, ...
                'Color', colors(i,:), 'LineWidth', 2)
        end
    end
    if j==bestChannels(c); textColor='red'; else; textColor='black'; end
    text(double(xcoords(j)*xSpacing+spkWindowInds(1)), ...
        ycoords(j)*ySpacing, ...
        num2str(j), 'Color', textColor)
end
set(gca, 'visible', 'off')
    
    
    