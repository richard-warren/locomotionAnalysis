% all units, all tracks on ccf

% settings
rasterize3d = false;  % this stretches the image a little atm


unitInfo = getUnitInfo('nucleiOnly', false);
ccf = loadCCF();
paper2_config;
load(fullfile(getenv('OBSDATADIR'), 'histology', '0_ephysHistoData', 'ephysHistoTable.mat'), ...
    'ephysHistoTable')

nucbins = ismember(unitInfo.nucleus, {'fastigial', 'interpositus', 'dentate'});

% get colors for scatter points
[~, cind] = ismember(unitInfo.nucleus, {'dentate', 'interpositus', 'fastigial'});
cind(cind==0) = 4;
scatColors = [cfg.nucleusColors; 0 0 0];
scatColors = scatColors(cind,:);
colorRep = repelem(cfg.nucleusColors,2,1);


% 3D
close all
figure('color', 'white', 'position', [231.00 529.00 319.00 275.00], 'menubar', 'none');

if rasterize3d; props = get(gca); end
plotLabels3D(ccf.coarseLabels==2, 'method', 'contours', 'colors', [0 0 0], 'smoothing', 5, 'contourAlpha', .05, ...
    'apGrid', ccf.ap, 'dvGrid', ccf.dv, 'mlGrid', ccf.ml, 'downSampling', 2, 'slices', {'dv', 'ap'}, 'nLines', 30);
plotLabels3D(ccf.labels, 'method', 'boundary', 'colors', colorRep, ...
    'surfArgs', {'edgecolor', 'none', 'facealpha', .4}, 'apGrid', ccf.ap, 'dvGrid', ccf.dv, 'mlGrid', ccf.ml);


% add probe tracks
for i = 1:height(ephysHistoTable)  % loop over shanks for all recordings
    
    mouse = ephysHistoTable.mouseID{i};
    
    % get entry point of shank in original histo coordinates
    surfaceLoc = ephysHistoTable.BS_crossPoints(i,:);
    
    % get location of lowest unit
    unitLocs = ephysHistoTable.GC_Points{i};
    
    if ~isempty(unitLocs)
        [~, maxInd] = max(unitLocs(:,2));
        unitLoc = unitLocs(maxInd,:);

        % combine points for line
        pts = [surfaceLoc; unitLoc];

        % load linear transformation from histo to ccf coordinates for this side of brain
        load(fullfile(getenv('SSD'), 'paper2', 'histo', 'registration', ...
            [mouse '_registration.mat']), 'tforms');
        load(fullfile(getenv('SSD'), 'paper2', 'histo', 'histoLabels', [mouse '_histoLabels.mat']), 'scaling', 'ap');
        tform = tforms.(ephysHistoTable.side{i});

        % convert pts to pixels (tform is from histo pixels to ccf pixels)
        pts = pts / 1000;
        pts(:,[1 3]) = pts(:,[1 3]) * 500 * scaling;
        pts(:,2) = pts(:,2) / diff(ap(1:2));

        % apply transformation!
        pts = cat(2, pts, [1 1]') * tform;
        pts = pts(:, 1:3) * .025;  % from ccf pixels to ccf mm

        % aaaaaaaaaand plot
        plot3(pts(:,1), pts(:,2), pts(:,3), 'Color', [0 0 0 .3], 'LineWidth', 1)
    end
    
end

% scatter units
scatter3(unitInfo.ccfMm(nucbins,1), unitInfo.ccfMm(nucbins,2), unitInfo.ccfMm(nucbins,3), ...
        5, 'black', 'filled', 'MarkerFaceAlpha', 1);
scatter3(unitInfo.ccfMm(~nucbins,1), unitInfo.ccfMm(~nucbins,2), unitInfo.ccfMm(~nucbins,3), ...
        5, 'black', 'MarkerEdgeAlpha', .6);


% add axis labels
len = .5;  % mm
args = {'LineWidth', 2, 'color', get(gca, 'XColor')};
axloc = [4 11 7];  % (mm) ml ap dv
plot3(axloc(1)+[0 len], axloc(2)+[0 0], axloc(3)+[0 0], args{:});  text(axloc(1)+len*2, axloc(2), axloc(3), 'ML');
plot3(axloc(1)+[0 0], axloc(2)+[0 -len], axloc(3)+[0 0], args{:}); text(axloc(1), axloc(2)-len*2, axloc(3), 'AP');
plot3(axloc(1)+[0 0], axloc(2)+[0 0], axloc(3)+[0 -len], args{:}); text(axloc(1), axloc(2), axloc(3)-len*2, 'DV');

% fancify
set(gca, 'visible', 'off', 'View', [45 30])

if rasterize3d
    frame = getframe(gca);
    frame = frame2im(frame);
    cla
    for p = fieldnames(props)'; try; set(gca, p{1}, props.(p{1})); catch; end; end
    image(frame)
    set(gca, 'visible', 'off')
end

% save
set(gcf, 'Renderer', 'painters')  % ensures export is not rasterized
saveas(gcf, 'E:\lab_files\paper2\paper_figures\matlab\ccf_allunits3D.svg')

%% 2D projections
plotUnitsOnCcf(unitInfo, 'colors', scatColors, 'figpos', [767.00 532.00 253.00 400.00])
saveas(gcf, 'E:\lab_files\paper2\paper_figures\matlab\ccf_allunits2D.svg')





























