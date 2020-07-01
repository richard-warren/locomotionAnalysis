% function aggregatePlots()

s.predictors = {'reward_normal'};
s.binNum = 100;
s.eventLims = [-.5 1];
s.epochLims = [-.5 1];
s.contLims = [];  % determine this either: manually // take median of percentile lims across mice...

% todo: determine total number of cells in advance? // figure out xLims in
% advance or manually define per predictor... // each heatmap has table
% with metadata, eg MI, correlation, std, x confidence...


% find sessions
ephysInfo = readtable(fullfile(getenv('OBSDATADIR'), 'spreadSheets', 'ephysInfo.xlsx'));
sessions = ephysInfo.session(ephysInfo.include==1);
clear ephysInfo

% initializations
xEvent = linspace(s.eventLims(1), s.eventLims(2), s.binNum);
xEpoch = linspace(s.epochLims(1), s.epochLims(2), s.binNum);
xCont = [];
for i = 1:length(s.predictors); heatmaps.(s.predictors{i}) = nan(0,s.binNum); end  % a sad hack


% for all sessions
for i = 1:length(sessions)
    try
        fprintf('%s: adding responses from session...\n', sessions{i})
        neuralData = load(fullfile(getenv('OBSDATADIR'), 'sessions', sessions{i}, 'neuralData.mat'));
        load(fullfile(getenv('OBSDATADIR'), 'sessions', sessions{i}, 'modelling', 'responses.mat'));

        % for all cells
        for j = 1:length(responses)
            cellResponses = responses(j).responses;

            % for each predictor
            for k = 1:length(s.predictors)

                if cellResponses{s.predictors{k}, 'type'}=='event'
                    response = cellResponses{s.predictors{k}, 'response'}{1};
                    response = (response-nanmean(neuralData.spkRates(j,:))) / nanstd(neuralData.spkRates(j,:));  % z score

                    xOriginal = linspace(cellResponses{s.predictors{k}, 'xLims'}(1), ...
                                         cellResponses{s.predictors{k}, 'xLims'}(2), length(response));
                    heatmaps.(s.predictors{k})(end+1,:) = interp1(xOriginal, nanmean(response,1), xEvent);
                end
            end
        end
    catch
        disp('  skiping...')
    end
end

%% test plot

close all;
figure('color', 'white', 'position', [1415.00 156.00 361.00 764.00]); hold on
colormap hot
imagesc(xEvent, 1:size(heatmaps.reward_normal,1), heatmaps.reward_normal)
plot([0 0], [1, size(heatmaps.reward_normal,1)], 'color', 'blue')






