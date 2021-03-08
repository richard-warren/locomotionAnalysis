%% all units, all tracks on ccf

data = getUnitInfo(false);
ccf = loadCCF();
paper2_config;
load(fullfile(getenv('OBSDATADIR'), 'histology', '0_ephysHistoData', 'ephysHistoTable.mat'), ...
    'ephysHistoTable')

% get colors for scatter points
[~, cind] = ismember(data.nucleus, {'interpositus', 'dentate', 'fastigial'});
cind(cind==0) = 4;
scatColors = [cfg.nucleusColors; 0 0 0];
scatColors = scatColors(cind,:);

%%



% 3D
close all
figure('color', 'white', 'position', [429.00 574.00 770.00 701.00]);

% wires, cerebellum only outline
plotLabels3D(ccf.coarseLabels==2, 'method', 'contours', 'colors', [0 0 0], 'smoothing', 5, 'contourAlpha', .05, ...
    'apGrid', ccf.ap, 'dvGrid', ccf.dv, 'mlGrid', ccf.ml, 'downSampling', 2, 'slices', {'dv', 'ap'}, 'nLines', 30);

% % flat, cerebellum only outline
% plotLabels3D(ccf.coarseLabels==2, 'method', 'boundary', 'colors', [0 0 0], 'smoothing', 5, ...
%     'surfArgs', {'edgecolor', 'none', 'facealpha', .05}, 'apGrid', ccf.ap, 'dvGrid', ccf.dv, 'mlGrid', ccf.ml, 'downSampling', 2);

% % flat, entire brain outline
% plotLabels3D(ccf.coarseLabels, 'method', 'boundary', 'colors', zeros(3,3), 'smoothing', 5, ...
%     'surfArgs', {'edgecolor', 'none', 'facealpha', .05}, 'apGrid', ccf.ap, 'dvGrid', ccf.dv, 'mlGrid', ccf.ml, 'downSampling', 2);

plotLabels3D(ccf.labels, 'method', 'boundary', 'colors', repelem(cfg.nucleusColors,2,1), ...
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


% add units
scatter3(data.ccfMm(:,1), data.ccfMm(:,2), data.ccfMm(:,3), ...
        20, 'black', 'filled', 'MarkerFaceAlpha', 1);
set(gca, 'visible', 'off', 'View', [45 30])
set(gcf, 'Renderer', 'painters')  % ensures export is not rasterized
saveas(gcf, 'E:\lab_files\paper2\paper_figures\matlab_exports\ccf_allunits', 'svg')













%%