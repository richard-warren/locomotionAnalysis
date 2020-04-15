function sesPlotRick(data, varargin)

% given matrix where first dim is mouse number, and second dim is session
% number, plots variavle for all mice as function of session number, eg to
% watch learning over time

% settings
s.colors = 'hsv';
s.scatSize = 20;
s.alpha = .4;
s.meanColor = [.2 .2 .2];
s.xvals = [];
s.xlabel = [];
s.ylabel = [];
s.lineWidth = 4;
s.compareTo = [];  % run significance tests comparing all sessions to averaged session numbers in 'compareTo'


% initializations
if exist('varargin', 'var'); for i = 1:2:length(varargin); s.(varargin{i}) = varargin{i+1}; end; end  % reassign settings contained in varargin
if ischar(s.colors); s.colors= eval([s.colors '(size(data,1))']); end
if isempty(s.xvals); s.xvals = 1:size(data,2); end
hold on

% plot thin lines
for i = 1:size(data,1)
    plot(s.xvals, data(i,:), 'LineWidth', 1, 'Color', [s.colors(i,:) s.alpha])
    scatter(s.xvals, data(i,:), s.scatSize, s.colors(i,:), 'filled', 'MarkerFaceAlpha', s.alpha)
end

% plot mean and std
plot(s.xvals, nanmean(data,1), 'LineWidth', s.lineWidth, 'Color', s.meanColor);
for i = 1:size(data,2)
    std = nanstd(data(:,i));
    line([s.xvals(i) s.xvals(i)], [-std std]+nanmean(data(:,i)), 'Color', s.meanColor, 'lineWidth', s.lineWidth/2)
end

% pimp fig
if ~isempty(s.xlabel); xlabel(s.xlabel); end
if ~isempty(s.ylabel); ylabel(s.ylabel); end
set(gca, 'TickDir', 'out')

pause(.001)


if ~isempty(s.compareTo)
    bl = mean(data(:, s.compareTo),2);
    for i = s.compareTo(end)+1:size(data,2)
%         [~, p] = ttest(bl, data(:,i));
        [p] = signrank(bl, data(:,i));
        fprintf('session %i, p: %.5f\n', i, p);
    end
    fprintf('\n');
end



