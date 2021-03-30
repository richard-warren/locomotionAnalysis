function errorplot(x, mat, varargin)
% rick version of shadedErrorBars

s.error = 'std';  % 'std' or 'sem'
s.color = lines(1);
s.linewidth = 2;
s.alpha = .1;


% reassign settings passed in varargin
if exist('varargin', 'var'); for i = 1:2:length(varargin); s.(varargin{i}) = varargin{i+1}; end; end

% get mean and error bars
mn = nanmean(mat, 1);
if strcmp(s.error, 'std')
    st = nanstd(mat, 1);
elseif strcmp(s.error, 'sem')
    st = nanstd(mat, 1) / sqrt(size(mat,1));
end

% plot
plot(x, mn, 'LineWidth', s.linewidth, 'color', s.color);
patch([x fliplr(x)], [mn+st fliplr(mn-st)], s.color, ...
    'EdgeColor', 'none', 'FaceAlpha', s.alpha)

