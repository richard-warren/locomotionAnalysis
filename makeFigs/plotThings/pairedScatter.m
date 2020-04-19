function pairedScatter(data, varargin)

% create simple scatter plot for paired data along with paired significance
% test...

% settings
s.ylabel = '';
s.colors = [.15 .15 .15; .15 .15 .15];  % colors associated with each axis
s.axisProps = {};  % properties that get applied to axix
s.scatterSize = 20;
s.scatterColors = [.15 .15 .15];
s.lims = [];  % use same limits for x and y axis

% stats
s.test = 'ttest';                 % 'ttest' or 'signrank'
s.pThresh = [.05 .01 .001];       % !!! CURRENTLY MUST BE ORDERED FROM LARGEST TO SMALLEST!
s.symbols = {'*', '**', '***'};   % symbols associated with the pThresh values above (they will appear above the brackets connecting the conditions to be compared)
s.notSigText = 'n.s.';            % text to appear above brackets for not groups that do not significanctly differ


% reassign settings passed in varargin
if exist('varargin', 'var'); for i = 1:2:length(varargin); s.(varargin{i}) = varargin{i+1}; end; end
if ischar(s.scatterColors); s.scatterColors = eval([s.scatterColors '(dataDims(end))']); end  % set scatter colors if color is specified as a string


% scatter
scatter(data(1,:), data(2,:), s.scatterSize, s.scatterColors, 'filled'); hold on

% get axis limits
if isempty(s.lims)
    lims = [get(gca, 'xlim'); get(gca, 'ylim')];
    s.lims = [min(lims(:,1)) max(lims(:,2))];
end

% add unity line
unity = plot(s.lims, s.lims, 'LineWidth', 1, 'Color', [.8 .8 .8]);
uistack(unity, 'bottom')

% set axis properties
% set(gca, 'xlim', s.lims, 'ylim', s.lims, 'XColor', s.colors(1,:), 'YColor', s.colors(2,:), s.axisProps{:})
set(gca, 'xlim', s.lims, 'ylim', s.lims, s.axisProps{:})
daspect([1 1 1])

% add labels
xlabel(s.xlabel, 'Color', s.colors(1,:))
ylabel(s.ylabel, 'Color', s.colors(2,:))
title(s.title, 'FontWeight', 'normal')



% run statistics
switch s.test
    case 'ttest'
        [~, p] = ttest(data(1,:), data(2,:));
    case 'signrank'
        p = signrank(data(1,:), data(2,:));
end

% determine whether significance is reached
pInd = find(p<s.pThresh, 1, 'last');
isSig = ~isempty(pInd);

% set significance dependent parameters
if isSig
    t = s.symbols{pInd};
    c = s.colors(2,:);
else
    t = s.notSigText;
    c = [.15 .15 .15];
end

% print results to command line
fprintf('%s-> p = %.2d = %.5f %s\n', s.test, p, p, t)
fprintf('difference-> %.2f\n\n', nanmean(data(1,:)) - nanmean(data(2,:)))

% add text above bracket
text(s.lims(1) + range(s.lims)*.5, s.lims(2), t, 'HorizontalAlignment', 'center', 'Color', c, 'VerticalAlignment', 'top')


