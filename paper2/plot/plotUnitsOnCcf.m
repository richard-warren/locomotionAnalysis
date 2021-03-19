function plotUnitsOnCcf(unitInfo, varargin)

% plot units on ccf, with units colored by by row sin unitInfo // unitInfo
% created by getUnitInfo()

% temp
% unitInfo = getUnitInfo();

% settings
s.colors = [0 0 0];
s.scalebar = .25;
s.half = true;
s.views = {'ap', 'dv'};
s.figpos = [767.00 532.00 253.00 200*length(s.views)];

% inits
if exist('varargin', 'var'); for i = 1:2:length(varargin); s.(varargin{i}) = varargin{i+1}; end; end  % reassign settings passed in varargin
ccf = loadCCF();
paper2_config;
nucColors = repelem(cfg.nucleusColors,2,1);
close all
figure('color', 'white', 'position', s.figpos, 'menubar', 'none')

if s.half
    leftIds = find(contains(ccf.nuclei, 'Left'));
    labels = ccf.labels .* uint8(ismember(ccf.labels, leftIds));
    
    midline = (ccf.ml(end)-ccf.ml(1)) / 2;
    ml = unitInfo.ccfMm(:,1);
    bins = ml>midline;
    ml(bins) = -(ml(bins) - midline) + midline;
    unitInfo.ccfMm(:,1)= ml;
else
    labels = ccf.labels;
end

% 2D projections
dims = {'ml', 'ap', 'dv'};  % names of axes (eg names of columns in unitInfo.ccfMm)

for i = 1:length(s.views)
    subplot(length(s.views),1,i); hold on
    inds = find(~strcmp(dims, s.views{i}));  % indices into 3D dimenions for x and y axes of this plot
    
    plotLabels2D(labels, 'dim', s.views{i}, ...
        'colors', nucColors, 'apGrid', ccf.ap, 'dvGrid', ccf.dv, 'mlGrid', ccf.ml, ...
        'patchArgs', {'FaceColor', 'none'})
    
    scatter(unitInfo.ccfMm(:, inds(1)), unitInfo.ccfMm(:, inds(2)), ...
        15, s.colors, 'filled', 'MarkerFaceAlpha', .6);
    
    % adjust axis limits
    for j = 1:2  % for x, then y axis
        otherdims = find(~strcmp(dims, dims{inds(j)}));  % 3D dims other than the current axis
        axsub = ccf.(dims{inds(j)})(any(permute(labels, [3 1 2]), otherdims));  % axis values for portion of axis containing nuclei labels (permute accounts for the fact that order of dimensions in 'labels' is [ap, dv, ml])
        buffer = range(axsub)*.1;                        % amount to pad the axis by on either side of the labeled portion
        lims = [axsub(1)-buffer axsub(end)+buffer];
        if j==1; set(gca, 'xlim', lims); else; set(gca, 'ylim', lims); end
    end
    
    % scale bars
    xflip = strcmp(get(gca, 'XDir'), 'reverse');
    yflip = strcmp(get(gca, 'YDir'), 'reverse');
    
    xlims = xlim; ylims = ylim;
    xlabel = get(get(gca, 'XLabel'), 'string');
    ylabel = get(get(gca, 'YLabel'), 'string');
    
    x = xlims(1+xflip) + [0 0 s.scalebar*(1-2*xflip)];
    y = ylims(1+yflip) + [s.scalebar*(1-2*yflip) 0 0];
    plot(x, y, 'LineWidth', 2, 'color', get(gca, 'XColor'));
    text(x(1), mean(y(1:2)), ylabel, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle')
    text(mean(x(2:3)), y(2), xlabel, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top')
    set(gca, 'visible', 'off')
end


