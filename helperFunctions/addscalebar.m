function addscalebar(scale, xy, txt)
% adds scale bar to lower left corner of plot // xy is either 'x' or 'y'

xlims = xlim;
ylims = ylim;

if strcmp(xy, 'x')
    x = xlims(1) + [0 scale];
    y = ylims(1) + [0 0];
    
    if exist('txt', 'var')
%         keyboard
        text(mean(x), y(1), txt, ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'top')
    end
    
elseif strcmp(xy, 'y')
    x = xlims(1) + [0 0];
    y = ylims(1) + [0 scale];
    
    if exist('txt', 'var')
        text(mean(x), mean(y), txt, 'Rotation', 90, ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom')
    end
end

plot(x, y, 'LineWidth', 2, 'Color', 'black')

