function addscalebar(scale, xy, txt, font, fontsize, color)
% adds scale bar to lower left corner of plot // xy is either 'x' or 'y'

if nargin<4; font = 'Times'; end
if nargin<5; fontsize = 8; end
if nargin<6; color = get(gca, 'XColor'); end

xlims = xlim;
ylims = ylim;

if strcmp(xy, 'x')
    x = xlims(1) + [0 scale];
    y = ylims(1) + [0 0];
    
    if exist('txt', 'var')
%         keyboard
        text(mean(x), y(1), txt, 'FontName', font, 'FontSize', fontsize, ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', ...
            'Color', color)
    end
    
elseif strcmp(xy, 'y')
    x = xlims(1) + [0 0];
    y = ylims(1) + [0 scale];
    
    if exist('txt', 'var')
        text(mean(x), mean(y), txt, 'Rotation', 90, 'FontName', font, ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
            'FontSize', fontsize, 'Color', color)
    end
end

plot(x, y, 'LineWidth', 2, 'Color', color)

