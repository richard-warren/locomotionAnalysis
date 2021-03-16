function addylabel(label)
% add x label using 'text' function // useful when need to hide y axis

xlims = xlim;
ylims = ylim;
text(xlims(1) - .1*range(xlims), mean(ylims), label, 'Rotation', 90, ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom')  % ylabel