function registerBrain(mouse)

% todo: save table with transformation, unit_ids, locations in pixels, mm,
% ccf, histo, etc...

% --------
% REGISTER
% --------

% load ccf, mouse brain, and cell locations
fprintf('registering %s brain to allen brain common coordinate framework...\n', mouse)
ccf = loadCCF();
data = load(fullfile(getenv('SSD'), 'paper2', 'histo', 'histoLabels', [mouse '_histoLabels.mat']));
load(fullfile(getenv('OBSDATADIR'), 'histology', '0_ephysHistoData', 'ephysHistoTable.mat'))
cellLocations = ephysHistoTable.GC_Points(strcmp(ephysHistoTable.mouseID, mouse));
cellLocations = cat(1, cellLocations{:});  % ml ap dv


% find transformation from histo (pixels) to ccf (pixels)

% 1) histo to cropped histo coords
apOffset = find(any(data.labels, [2 3]), 1, 'first');
dvOffset = find(any(data.labels, [1 3]), 1, 'first');
mlOffset = find(any(data.labels, [1 2]), 1, 'first');
T1 = eye(4);
T1(end,1:3) = -[dvOffset apOffset mlOffset];  % dv ap ml

% 2) cropped histo to resized coords
apScale = sum(any(ccf.labels, [2 3])) / sum(any(data.labels, [2 3]));
dvScale = sum(any(ccf.labels, [1 3])) / sum(any(data.labels, [1 3]));
mlScale = sum(any(ccf.labels, [1 2])) / sum(any(data.labels, [1 2]));
T2 = diag([dvScale apScale mlScale 1]);  % dv ap ml

% 3) cropped histo to cropped ccf
[optimizer, metric] = imregconfig('monomodal');
optimizer.MaximumIterations = 200;
labelsCcfCropped = ccf.labels(any(ccf.labels, [2 3]), any(ccf.labels, [1 3]), any(ccf.labels, [1 2]));
labelsCropped = data.labels(any(data.labels, [2 3]), any(data.labels, [1 3]), any(data.labels, [1 2]));
labelsCropped = imresize3(labelsCropped, size(labelsCcfCropped), 'nearest');
tform = imregtform(labelsCropped, labelsCcfCropped, 'affine', optimizer, metric, 'DisplayOptimization', false);
T3 = tform.T;

% 4) cropped ccf to normal ccf
apOffset = find(any(ccf.labels, [2 3]), 1, 'first');
dvOffset = find(any(ccf.labels, [1 3]), 1, 'first');
mlOffset = find(any(ccf.labels, [1 2]), 1, 'first');
T4 = eye(4);
T4(end,1:3) = [dvOffset apOffset mlOffset];  % dv ap ml

% full transform
T = T1 * T2 * T3 * T4;
tform = affine3d(T);
warped = imwarp(data.labels, tform, 'OutputView', imref3d(size(ccf.labels)), 'interp', 'nearest');

% convert cell locations to pixels
cellLocationsHistoMm = cellLocations / 1000;
cellLocationsHistoPixels = cellLocationsHistoMm;  % ml ap dv
cellLocationsHistoPixels(:,[1 3]) = cellLocationsHistoPixels(:,[1 3]) * 500 * data.scaling;
cellLocationsHistoPixels(:,2) = cellLocationsHistoPixels(:,2) / diff(data.ap(1:2));

% apply transormation to cell locations
% (columns of T act on (dv, ap, ml), but locations are (ml, ap, dv), so
% first we create T_cells, which swaps the order of the diagonal entries
% and the bottom row)
T_cells = T(:,[3 2 1 4]);
T_cells(1:3,1:3) = flipud(T_cells(1:3,1:3));
cellLocationsCcfPixels = [cellLocationsHistoPixels, ones(size(cellLocations,1),1)] * T_cells;
cellLocationsCcfPixels = cellLocationsCcfPixels(:,1:3);
cellLocationsCcfMm = cellLocationsCcfPixels * .025;


% ----
% PLOT
% ----

% 3d plot
colors = repelem(lines(3),2,1);
figure('color', 'white', 'position', [79.00 48.00 1794.00 928.00]); hold on
plotLabels3D(ccf.labels, 'surfArgs', {'FaceAlpha', .1, 'edgealpha', .2}, 'colors', colors, ...
    'apGrid', ccf.ap, 'dvGrid', ccf.dv, 'mlGrid', ccf.ml);
plotLabels3D(warped, 'colors', colors, 'apGrid', ccf.ap, 'dvGrid', ccf.dv, 'mlGrid', ccf.ml);
scatter3(cellLocationsCcfMm(:,1), cellLocationsCcfMm(:,2), cellLocationsCcfMm(:,3), ...
    50, [0 0 0], 'filled');
savefig(fullfile(getenv('OBSDATADIR'), 'figures', 'brainRegistration', [mouse '_3D.fig']))

% 2d plots
figure('color', 'white', 'position', [79.00 48.00 1794.00 928.00]);
views = {'ap', 'dv', 'ml'};
inds = {[1 3], [1 2], [2 3]};  % cellLocations inds for 'ap' and 'ml' views

for i = 1:3
    % zoomed out
    subplot(3,3,1 + (i-1)*3); hold on
    title('registered')
    plotLabels2D(warped, 'dim', views{i}, ...                   % registered
        'colors', colors, 'apGrid', ccf.ap, 'dvGrid', ccf.dv, 'mlGrid', ccf.ml)  
    plotLabels2D(ccf.coarseLabels, 'dim', views{i}, ...         % cerebellum outline
        'patchArgs', {'facecolor', [0 0 0], 'facealpha', .05, 'edgecolor', [0 0 0], 'edgealpha', .4}, ...
        'apGrid', ccf.ap, 'dvGrid', ccf.dv, 'mlGrid', ccf.ml, 'smoothing', 10)
    scatter(cellLocationsCcfMm(:,inds{i}(1)), cellLocationsCcfMm(:,inds{i}(2)), ...
        30, 'black', 'filled', 'MarkerFaceAlpha', .6);

    % original
    subplot(3,3,2 + (i-1)*3); hold on
    title('original')
    plotLabels2D(data.labels, 'dim', views{i}, ...
        'colors', colors, 'apGrid', data.ap, 'dvGrid', data.dv, 'mlGrid', data.ml)
    scatter(cellLocationsHistoMm(:,inds{i}(1)), cellLocationsHistoMm(:,inds{i}(2)), ...
        30, 'black', 'filled', 'MarkerFaceAlpha', .6);
    
    % registered
    subplot(3,3,3 + (i-1)*3); hold on
    title('registered')
    plotLabels2D(warped, 'dim', views{i}, ...                   % registered
        'colors', colors, 'apGrid', ccf.ap, 'dvGrid', ccf.dv, 'mlGrid', ccf.ml)  
    plotLabels2D(ccf.labels, 'dim', views{i}, ...               % ccf
        'patchArgs', {'facecolor', 'none', 'edgecolor', [0 0 0], 'edgealpha', .4}, ...
        'apGrid', ccf.ap, 'dvGrid', ccf.dv, 'mlGrid', ccf.ml)
    scatter(cellLocationsCcfMm(:,inds{i}(1)), cellLocationsCcfMm(:,inds{i}(2)), ...
        30, 'black', 'filled', 'MarkerFaceAlpha', .6);
end

saveas(gcf, fullfile(getenv('OBSDATADIR'), 'figures', 'brainRegistration', [mouse '_2D.png']))
save(fullfile(getenv('SSD'), 'paper2', 'histo', 'registration', [mouse '_registration.mat']), ...
    'cellLocationsCcfMm', 'cellLocationsCcfPixels')
disp('all done!')
























