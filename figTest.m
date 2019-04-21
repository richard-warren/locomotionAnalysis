
close all; figure; plot(1:100, rand(5,100))
set(gca, 'ylim', [-1 2])

file = fullfile(getenv('OBSDATADIR'), 'papers', 'paper1', 'figures', 'matlabFigs', ...
        'test');
saveas(gcf, file, 'epsc');
savefig(gcf, file);