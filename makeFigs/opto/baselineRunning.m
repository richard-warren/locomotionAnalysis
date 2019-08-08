% plots the results of experiments in which the effects of
% light on baseline running (with no obstacles) are tested
% !!! modify to automatically accept powers used in exp...


% settings
% sessions = {'190805_000', '190805_001', '190805_002', '190806_000', '190806_001', '190806_002'}; % old sessions with short light duration and old stim intensities
% powers = [0, .25, .5, 1];  % laser percentages used in experiment
sessions = {'190807_000', '190807_001', '190807_002'};
powers = [0, 0.12, 0.35, 1.0];
normalizeVel = true;  % whether to subtract vel at moment of obs on

velWindow = .1;
velTimes = [-.5 2];
wheelFs = 1000;
errorFcn = @(x) nanstd(x) /sqrt(size(x,1));
rampUpTime = .2;
sinFreq = 40;
sinHgt = .15;
yLims = [-.8 .8];
yLimsNormalized = [-.6 .2];  % yLims to use when velocity is normalized


% initializations
if normalizeVel; yLims = yLimsNormalized; end
sessionInfo = readtable(fullfile(getenv('OBSDATADIR'), 'spreadSheets', 'experimentMetadata.xlsx'), 'Sheet', 'optoNotes');
sessionInfo = sessionInfo(ismember(sessionInfo.session, sessions),:); % remove empty rows, not included sessions
mice = unique(sessionInfo.mouse);
brainRegions = unique(sessionInfo.brainRegion);
data = struct();
rowInd = 1;
x = velTimes(1):(1/wheelFs):velTimes(2);
colors = [0 0 0; winter(length(powers)-1)]; % scale linearly from black to color
lightOnInd = find(x==0);

% collect data all data in big ol table
for i = 1:length(sessions)
    disp(sessions{i})
    
    % load session data
    load(fullfile(getenv('OBSDATADIR'), 'sessions', sessions{i}, 'runAnalyzed.mat'), ...
        'wheelPositions', 'wheelTimes', 'obsOnTimes', 'obsOffTimes', 'targetFs')
    load(fullfile(getenv('OBSDATADIR'), 'sessions', sessions{i}, 'run.mat'), 'stimulus')
    wheelVel = getVelocity(wheelPositions, velWindow, targetFs);
    stimulusTimes = stimulus.times;
    stimulus = stimulus.values;
    sessionBin = strcmp(sessionInfo.session, sessions{i});
    
    % loop over trials
    for j = 1:length(obsOnTimes)
        
        % get trial light condition
        trialStim = stimulus(stimulusTimes>obsOnTimes(j) & stimulusTimes<obsOffTimes(j));
        trialPower = max(trialStim)/5;  % peak signal is trialPower fraction of 5V max
        [minDif, minInd] = min(abs(powers-trialPower));  % find closest power in powers, defined above
        if minDif<.01; trialPower = powers(minInd); else; trialPower = nan; end  % trialPower is nan if close value is not in powers, defined above
        
        % get vel vs time
        startInd = find(wheelTimes>=obsOnTimes(j)+velTimes(1), 1, 'first');
        if length(wheelVel)>=startInd+length(x)-1
            trialVel = wheelVel(startInd:startInd+length(x)-1);
        else
            trialVel = nan(1,length(x)); 
        end
        
        % store in data struct
        data(rowInd).session = sessions{i};
        data(rowInd).mouse = sessionInfo.mouse{sessionBin};
        data(rowInd).region = sessionInfo.brainRegion{sessionBin};
        data(rowInd).trial = j;
        data(rowInd).lightPower = trialPower;
        data(rowInd).trialVel = trialVel;
        data(rowInd).obsOnDurations = obsOffTimes(j) - obsOnTimes(j);
        rowInd = rowInd + 1;
        
    end
end


% average across each mouse
dataMatrix = nan(length(brainRegions), length(powers), length(mice), length(x));
for i = 1:length(brainRegions)
    for j = 1:length(powers)
        for k = 1:length(mice)            
            bins = strcmp({data.region}, brainRegions{i}) & ...
                   [data.lightPower]==powers(j) & ...
                   strcmp({data.mouse}, mice{k});
            dataSub = cat(1, data(bins).trialVel);
            if normalizeVel; dataSub = dataSub - dataSub(:, lightOnInd); end
            dataMatrix(i,j,k,:) = nanmean(dataSub, 1); % get avg for this brain region, for this light power, for this mouse
        end
    end
end


%% plot means
close all;
figure('color', 'white', 'menubar', 'none', 'Position', [1960 156 730 757]);


for i = 1:length(brainRegions)
    
    subplot(length(brainRegions),1,i); hold on
    title(brainRegions{i})
    
    for j = 1:length(powers)             
        lineData = squeeze(dataMatrix(i,j,:,:));
        shadedErrorBar(x, lineData, {@nanmean, errorFcn}, ...
            'lineprops', {'linewidth', 5, 'color', colors(j,:)}, 'patchSaturation', .1);
%         plot(x, nanmean(lineData, 1), 'linewidth', 5, 'color', colors(j,:));
    end
    
    % plot sin wave
    sinX = x(x>=0);
    sinWave = (sin(sinX*2*pi*sinFreq)+1) * sinHgt/2;
    rampBins = sinX<=rampUpTime;
    sinWave(rampBins) = sinWave(rampBins) .* linspace(0,1,sum(rampBins));
    rampDownStart = median([data.obsOnDurations]); % median time at which stim starts to ramp down
    rampBins = sinX>=rampDownStart & sinX<(rampDownStart+rampUpTime);
    sinWave(rampBins) = sinWave(rampBins) .* linspace(1,0,sum(rampBins));
    sinWave = sinWave(1:find(rampBins,1,'last'));
    sinX = sinX(1:find(rampBins,1,'last'));
    sinWave = sinWave + yLims(1);
    plot(sinX, sinWave, 'LineWidth', 2, 'Color', mean(colors,1));
    
    set(gca, 'yLim', yLims, 'TickDir', 'out');
end
xlabel('time from laser (s)')
if yLimsNormalized; yLabel='\Delta velocity (m/s)'; else; yLabel='velocity (m/s)'; end
ylabel(yLabel)


% plot per mouse
figure('color', 'white', 'menubar', 'none', 'Position', [1946 26 1032 443]);


for i = 1:length(brainRegions) 
    for m = 1:length(mice)
        plotInd = (i-1)*length(mice) + m;
        subplot(length(brainRegions), length(mice), plotInd); hold on
        title([brainRegions{i} ': ' mice{m}])
        
        line([0 0], yLims)
        for j = 1:length(powers)
            lineData = squeeze(dataMatrix(i,j,m,:));
            plot(x, lineData, 'linewidth', 3, 'color', colors(j,:));
        end
        
        set(gca, 'yLim', yLims, 'TickDir', 'out');
    end
end














