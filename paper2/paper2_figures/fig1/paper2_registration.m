mouse = 'cer28';
slice = 10;
paper2_config;

d = load(['E:\lab_files\paper2\histo\histoLabels\' mouse '_histoLabels.mat']);
img = imread(['E:\lab_files\paper2\histo\scans\' mouse '_slice' num2str(slice) '.png']);
img = imresize(img, d.scaling);
c = repelem(cfg.nucleusColors, 2, 1);
load(fullfile(getenv('OBSDATADIR'), 'histology', '0_ephysHistoData', 'ephysHistoTable.mat'), 'ephysHistoTable')
probes = ephysHistoTable(strcmp(ephysHistoTable.mouseID, mouse), :);
clear ephysHistoTable

%% get probe lines

probes.locHisto = nan(height(probes), 2, 3);  % for each probe record xyz location of (1) entry point and (2) location of deepest unit
probes.locCcf = nan(height(probes), 2, 3);

% add probe tracks
for i = 1:height(probes)  % loop over shanks for all recordings
    
    
    surfaceLoc = probes.BS_crossPoints(i,:);  % get entry point of shank in original histo coordinates
    unitLocs = probes.GC_Points{i};  % get location of lowest unit
    
    if ~isempty(unitLocs)
        [~, maxInd] = max(unitLocs(:,2));
        unitLoc = unitLocs(maxInd,:);

        % combine points for line
        pts = [surfaceLoc; unitLoc] / 1000;
        probes.locHisto(i,:,:) = pts;

        % load linear transformation from histo to ccf coordinates for this side of brain
        load(fullfile(getenv('SSD'), 'paper2', 'histo', 'registration', ...
            [mouse '_registration.mat']), 'tforms');
        load(fullfile(getenv('SSD'), 'paper2', 'histo', 'histoLabels', [mouse '_histoLabels.mat']), 'scaling', 'ap');
        tform = tforms.(probes.side{i});

        % convert pts to pixels (tform is from histo pixels to ccf pixels)
%         pts = pts / 1000;
        pts(:,[1 3]) = pts(:,[1 3]) * 500 * scaling;
        pts(:,2) = pts(:,2) / diff(ap(1:2));

        % apply transformation!
        pts = cat(2, pts, [1 1]') * tform;
        pts = pts(:, 1:3) * .025;  % from ccf pixels to ccf mm
        probes.locCcf(i,:,:) = pts;
    end
end



%% histology with nucleus outlines
close all;
figure('color', 'black', 'position', [2.00 722.00 1278.00 634.00], 'menubar', 'none'); hold on

% zoomed out
subplot(1,2,1); hold on
image(d.ml, d.dv, img)
plotLabels2D(d.labels(slice,:,:), 'dvGrid', d.dv, 'mlGrid', d.ml, ...
    'colors', c, 'patchArgs', {'FaceColor', 'none', 'LineWidth', 3, 'LineStyle', ':'});
set(gca, 'Visible', 'off', 'color', 'black')
xlims = xlim;
plot(xlims(2)+[0 -1]-range(xlims)*.05, [5 5], 'Color', 'white', 'linewidth', 3)  % add scale bar
temp = getframe(gca); temp = frame2im(temp); imwrite(temp, 'E:\lab_files\paper2\paper_figures\matlab\histo_zoomout.png')



% zoomed in
subplot(1,2,2); hold on
xlims = [1.2 3.4]; ylims = [1.6 3.6];
image(d.ml, d.dv, img)
plotLabels2D(d.labels(slice,:,:), 'dvGrid', d.dv, 'mlGrid', d.ml, ..., ...
    'colors', c, 'patchArgs', {'FaceColor', 'none', 'LineWidth', 6, 'LineStyle', ':'});
set(gca, 'xlim', xlims, 'ylim', ylims, 'color', 'black')
set(gca, 'Visible', 'off')
plot(xlims(2)+[0 -.5]-range(xlims)*.02, ylims(2)+[0 0]-range(ylims)*.02, 'Color', 'white', 'linewidth', 5)  % plot scale bar
temp = getframe(gca); temp = frame2im(temp); imwrite(temp, 'E:\lab_files\paper2\paper_figures\matlab\histo_zoomin.png')


%% 3D view of all traced stuff

close all
figure('color', 'white', 'position', [483.00 331.00 797.00 1025.00], 'menubar', 'none'); hold on
transparency = .2;

x = probes.locHisto(:,:,1)'; y = probes.locHisto(:,:,3)';  z = probes.locHisto(:,:,2)'; % get probe locations

[X,Z] = meshgrid(d.ml, d.dv);

for i = 1:size(d.labels,1)
    
    Y = ones(size(X)) * d.ap(i);
    
    % plot purkinje cell layers and brain surface with surface()
    cmat = ones(size(d.labels,2), size(d.labels,3), 3);
    
    % the following lines incorporate the nuclei into surface code, which
    % is faster but doesn't allow the flexibility of patch()
%     for j = 1:6
%         bins = repmat(squeeze(d.labels(i,:,:))==j, 1, 1, 3);
%         cmat(bins) = repelem(c(j,:), 1, sum(bins(:))/3);
%     end
    
    cmat(repmat(squeeze(d.pcLayers(i,:,:) | d.surface(i,:,:)), 1, 1, 3)) = .5;
    transp = double(~all(cmat==1, 3)) * transparency;
    transp(1,1) = 1;  % hack // matlab appears to scale the alpha by the maximum value, so this ensures max(alpha)=1
    surface(X, Y, Z, cmat, 'edgecolor', 'none', 'cdatamapping', 'direct', ...
        'facecolor', 'flat', 'AlphaData', transp, 'FaceAlpha', 'flat');
    
    % use patch() to plot nuclei
    for j = 1:6
        mask = squeeze(d.labels(i,:,:)==j);
        if any(mask(:))
            [i1, i2] = ind2sub(size(mask), find(mask));
            dv = d.dv(i1); ml = d.ml(i2);
            k = boundary(dv', ml');
            dv = dv(k)'; ml = ml(k)'; ap = repelem(d.ap(i), length(dv));
            patch(ml, ap, dv, c(j,:), 'FaceAlpha', .1, 'EdgeColor', c(j,:), ...
                'LineWidth', 1, 'EdgeAlpha', .5)
        end
    end
end

plot3(x, z, y, 'color', [cfg.diiColor], 'LineWidth', 2)  % add probe traces

daspect([1 1 1])
set(gca, 'ZDir', 'reverse')
view(45, 30)

% add axis labels
len = .25;  % mm
args = {'LineWidth', 2, 'color', get(gca, 'XColor')};
axloc = [2 1 5];  % (mm) ml ap dv
plot3(axloc(1)+[0 len], axloc(2)+[0 0], axloc(3)+[0 0], args{:}); % text(axloc(1)+len*2, axloc(2), axloc(3), 'ML');
plot3(axloc(1)+[0 0], axloc(2)+[0 -len], axloc(3)+[0 0], args{:}); % text(axloc(1), axloc(2)-len*2, axloc(3), 'AP');
plot3(axloc(1)+[0 0], axloc(2)+[0 0], axloc(3)+[0 -len], args{:}); % text(axloc(1), axloc(2), axloc(3)-len*2, 'DV');

set(gca, 'visible', 'off')
saveas(gca, 'E:\lab_files\paper2\paper_figures\matlab\histo_tracing.png')













