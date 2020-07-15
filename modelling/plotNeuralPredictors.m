function plotNeuralPredictors(session, varargin)


% settings
s.dz = 6;  % (standard deviation) vertical separation
s.xWidth = 40;  % (seconds) range of x axis
s.xLims = [];
s.ommit = {'paw1LH_stride', 'paw2LF_stride', 'paw3RF_stride', 'paw4RH_stride', ...
    'reward', 'reward_surprise', 'reward_omission', 'reward_normal'};
s.predictors = [];  % can pass in predicotrs // otherwise it is loaded from disk
s.predictorList = {};


% initialiations
if exist('varargin', 'var'); for i = 1:2:length(varargin); s.(varargin{i}) = varargin{i+1}; end; end % reassign settings passed in varargin
if isempty(s.predictors)
    load(fullfile(getenv('OBSDATADIR'), 'sessions', session, 'modelling', 'predictors.mat'), 'predictors');
    s.predictors = predictors;
    clear predictors;
end
set(0, 'DefaultAxesTickLabelInterpreter', 'none')
figure('name', session, 'color', 'white', 'position', [125.00 93.00 1666.00 886.00]); hold on


% determine which predictors to show
if isempty(s.predictorList)
    allInds = find(~ismember(s.predictors.Properties.RowNames, s.ommit));
else
    [~, allInds] = ismember(s.predictorList, s.predictors.Properties.RowNames);  % if the user specifies the order of the predictors
    allInds = allInds(allInds>0);
end
s.predictors = s.predictors(allInds,:);  % restrict to desired predictors and shuffle to the correct order


% determine x axis
t = s.predictors{find(s.predictors.type=='continuous',1,'first'), 't'}{1};  % (assumes same time grid for all predictors!)
tMax = t(end);  % last time in first continuous predictor
if isempty(s.xLims); s.xLims = [0 s.xWidth] + randi(round(tMax-s.xWidth*2)); end

colors = lines(length(allInds));
y = 0 : s.dz : s.dz*height(s.predictors)-1;

% events
eventInds = find(s.predictors.type=='event');
plot(s.xLims, repmat(y(eventInds),2,1)', 'color', [.4 .4 .4])  % horizontal lines for each event predictor
for i = 1:length(eventInds)
    x = s.predictors{eventInds(i),'data'}{1};
    x = x(x>s.xLims(1) & x<s.xLims(2));
    scatter(x, repelem(y(eventInds(i)), length(x)), 20, colors(eventInds(i),:), 'filled')
end

% epoch
epochInds = find(s.predictors.type=='epoch');
plot(s.xLims, repmat(y(epochInds),2,1)', 'color', [.4 .4 .4])  % horizontal lines for each epoch predictor
for i = 1:length(epochInds)
    epoch = s.predictors{epochInds(i),'data'}{1};
    epoch = epoch(any(epoch>s.xLims(1) & epoch<s.xLims(2),2),:);  % only plot epochs within xLims
    if ~isempty(epoch)
        plot(epoch', repelem(y(epochInds(i)), 2), 'LineWidth', 5, 'color', colors(epochInds(i),:))
    end
end

% continuous
contBins = s.predictors.type=='continuous';
bins = t>s.xLims(1) & t<s.xLims(2);
allCont = cat(1, s.predictors{s.predictors.type=='continuous','data'}{:});
plot(t(bins)', zscore(allCont(:,bins)') + y(contBins), 'LineWidth', 1.5);



% fancify
set(gca, 'xlim', s.xLims, 'ytick', y, ...
    'YTickLabel', s.predictors.Properties.RowNames, ...
    'YLim', [y(1)-+s.dz y(end)+s.dz], 'TickDir', 'out', 'ydir', 'reverse')
end

