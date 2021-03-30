function addylabel(label, font, fontsize)
% add x label using 'text' function // useful when need to hide y axis

if nargin<2; font = 'Times'; end
if nargin<3; fontsize = 8; end

xlims = xlim;
ylims = ylim;
text(xlims(1) - .1*range(xlims), mean(ylims), label, 'Rotation', 90, ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
    'FontName', font, 'FontSize', fontsize)