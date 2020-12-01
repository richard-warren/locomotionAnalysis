function registerBrain(mouse)

% todo: save table with transformation, unit_ids, locations in pixels, mm,
% ccf, histo, etc...

% --------
% REGISTER
% --------

% settings
maxDistance = 80;  % (microns) cells within maxDistance of a nucleus are assigned to that nucleus



% load ccf, mouse brain, and cell locations
fprintf('registering %s to allen brain common coordinate framework... ', mouse)
ccf = loadCCF();
data = load(fullfile(getenv('SSD'), 'paper2', 'histo', 'histoLabels', [mouse '_histoLabels.mat']));
load(fullfile(getenv('OBSDATADIR'), 'histology', '0_ephysHistoData', 'ephysHistoTable.mat'), 'ephysHistoTable')

% make table for data storage
bins = strcmp(ephysHistoTable.mouseID, mouse);
cellLocations = ephysHistoTable.GC_Points(bins);
unitsPerSession = cellfun(@(x) size(x,1), cellLocations);
sessions = repelem(ephysHistoTable.session(bins), unitsPerSession);
unit_ids = ephysHistoTable.GC_ids(bins);
unit_ids = cat(1, unit_ids{:});
shank = repelem(ephysHistoTable.shankNum(bins), unitsPerSession);
nunits = sum(unitsPerSession);

% prepare nuclei names
for i = 1:length(data.nuclei); data.nuclei{i} = lower(erase(data.nuclei{i}, {'Left', 'Right'})); end  % remove Left, Right, and make lowercase


registration = table(sessions, unit_ids, shank, nan(nunits,3), nan(nunits,3), cell(nunits,1), ...
    'VariableNames', {'session', 'unit', 'shank', 'ccfMm', 'ccfPix', 'nucleus'});
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


% determine nucleus for each unit
% todo: determine if above nuclei, and therefore in cerebellar cortex

% first 'dilate' the nuclei to include cells that are really close
res = .025;  % resolution of resampled labels
dvgrid = data.dv(1):res:data.dv(end);
apgrid = data.ap(1):res:data.ap(end);
mlgrid = data.ml(1):res:data.ml(end);

maxPixelDif = round((maxDistance/1000) / res);  % units must be within this many pixels of a nucleus
labels = interp3(data.dv, data.ap, data.ml, double(data.labels), ...
    dvgrid', apgrid, mlgrid, 'nearest');  % interpolate such that x, y, and z resolution match and aren't too high
labelsDilated = imdilate(labels, strel('sphere', maxPixelDif));  % dilate
halo = labels==0 & labelsDilated>0;  % mask for pixels surrounding nuclei that are within maxPixelDif
labels(halo) = labelsDilated(halo);  % only use dilated labeld within the mask, becuse the dilation will cause weird things to happen within the mask (due to max operation across class labels)
labelsDilated = labels;

% uncomment to visualize 'expansion' of nuclei via dilation
% close all; figure('color', 'white', 'position', [2.00 722.00 1278.00 634.00]);
% plotLabels3D(data.labels, 'downSampling', 2, 'surfArgs', {'FaceAlpha', 1}, ...
%     'apGrid', data.ap, 'dvGrid', data.dv, 'mlGrid', data.ml); hold on;
% plotLabels3D(labels, 'downSampling', 2, 'surfArgs', {'FaceAlpha', .1}, ...
%     'apGrid', apgrid, 'dvGrid', dvgrid, 'mlGrid', mlgrid)

% ... then assign units to nuclei!
for i = 1:nunits
    mlInd = knnsearch(mlgrid', cellLocationsHistoMm(i,1));
    apInd = knnsearch(apgrid', cellLocationsHistoMm(i,2));
    dvInd = knnsearch(dvgrid', cellLocationsHistoMm(i,3));
    labelInd = labelsDilated(apInd, dvInd, mlInd);
    if labelInd>0
        registration{i, 'nucleus'} = {data.nuclei{labelInd}};
    else
        registration{i, 'nucleus'} = {'none'};
    end
end


% save registration
registration.ccfMm = cellLocationsCcfMm;
registration.ccfPix = cellLocationsCcfPixels;
save(fullfile(getenv('SSD'), 'paper2', 'histo', 'registration', [mouse '_registration.mat']), ...
    'mouse', 'registration')



% ----
% PLOT
% ----


% determine unit colors based on nuclei they are in
colorsTemp = lines(3);
scatColors = ones(nunits, 3);
for i = 1:3
    bins = strcmp(registration.nucleus, data.nuclei{i*2});  % a hack that assumes data.nuclei are format as {nuc1, nuc1, nuc2, nuc2, nuc3, nuc3}
    if any(bins); scatColors(bins, :) = repmat(colorsTemp(i,:), sum(bins), 1); end
end


% 3d plot
colors = repelem(lines(3),2,1);
figure('name', mouse, 'color', 'white', 'position', [79.00 48.00 1794.00 928.00]); hold on
plotLabels3D(ccf.labels, 'surfArgs', {'FaceAlpha', .1, 'edgealpha', .2}, 'colors', colors, ...
    'apGrid', ccf.ap, 'dvGrid', ccf.dv, 'mlGrid', ccf.ml);
plotLabels3D(warped, 'colors', colors, 'apGrid', ccf.ap, 'dvGrid', ccf.dv, 'mlGrid', ccf.ml);
scatter3(cellLocationsCcfMm(:,1), cellLocationsCcfMm(:,2), cellLocationsCcfMm(:,3), ...
    50, scatColors, 'filled', 'MarkerEdgeColor', [0 0 0], 'LineWidth', 3);
savefig(fullfile(getenv('SSD'), 'paper2', 'histo', 'registration', [mouse '_3D.fig']))

% 2d plots
figure('name', mouse, 'color', 'white', 'position', [79.00 48.00 1794.00 928.00]);
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
        30, scatColors, 'filled', 'MarkerFaceAlpha', .6, 'MarkerEdgeColor', 'black', 'LineWidth', 1);

    % original
    subplot(3,3,2 + (i-1)*3); hold on
    title('original')
    plotLabels2D(data.labels, 'dim', views{i}, ...
        'colors', colors, 'apGrid', data.ap, 'dvGrid', data.dv, 'mlGrid', data.ml)
    scatter(cellLocationsHistoMm(:,inds{i}(1)), cellLocationsHistoMm(:,inds{i}(2)), ...
        30, scatColors, 'filled', 'MarkerFaceAlpha', .6, 'MarkerEdgeColor', 'black', 'LineWidth', 1);
    
    % registered
    subplot(3,3,3 + (i-1)*3); hold on
    title('registered')
    plotLabels2D(warped, 'dim', views{i}, ...                   % registered
        'colors', colors, 'apGrid', ccf.ap, 'dvGrid', ccf.dv, 'mlGrid', ccf.ml)  
    plotLabels2D(ccf.labels, 'dim', views{i}, ...               % ccf
        'patchArgs', {'facecolor', 'none', 'edgecolor', [0 0 0], 'edgealpha', .4}, ...
        'apGrid', ccf.ap, 'dvGrid', ccf.dv, 'mlGrid', ccf.ml)
    scatter(cellLocationsCcfMm(:,inds{i}(1)), cellLocationsCcfMm(:,inds{i}(2)), ...
        30, scatColors, 'filled', 'MarkerFaceAlpha', .6, 'MarkerEdgeColor', 'black', 'LineWidth', 1);
end

saveas(gcf, fullfile(getenv('SSD'), 'paper2', 'histo', 'registration', [mouse '_2D.png']))
fprintf('all done!\n')
























