function sesPlotRick(data, opts)

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

% reassign settings contained in opts
if exist('opts', 'var'); for i = 1:2:length(opts); s.(opts{i}) = opts{i+1}; end; end


% initializations
if ischar(s.colors); s.colors= eval([s.colors '(size(data,1))']); end
if isempty(s.xvals); s.xvals = 1:size(data,2); end
hold on

% plot thin lines
for i = 1:size(data,1)
    plot(s.xvals, data(i,:), 'LineWidth', 1, 'Color', [s.colors(i,:) s.alpha])
    scatter(s.xvals, data(i,:), s.scatSize, s.colors(i,:), 'filled', 'MarkerFaceAlpha', s.alpha)
end

% plot mean and std
plot(s.xvals, nanmean(data,1), 'LineWidth', 4, 'Color', s.meanColor);
for i = 1:size(data,2)
    std = nanstd(data(:,i));
    line([s.xvals(i) s.xvals(i)], [-std std]+nanmean(data(:,i)), 'Color', s.meanColor, 'lineWidth', 2)
end

% pimp fig
if ~isempty(s.xlabel); xlabel(s.xlabel); end
if ~isempty(s.ylabel); ylabel(s.ylabel); end