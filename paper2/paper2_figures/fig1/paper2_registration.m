% inits

mouse = 'cer28';
slice = 10;       % slice to show for histology fig
paper2_config;

d = load(['E:\lab_files\paper2\histo\histoLabels\' mouse '_histoLabels.mat']);
img = imread(['E:\lab_files\paper2\histo\scans\' mouse '_slice' num2str(slice) '.png']);
img = imresize(img, d.scaling);
c = repelem(cfg.nucleusColors, 2, 1);
load(fullfile(getenv('OBSDATADIR'), 'histology', '0_ephysHistoData', 'ephysHistoTable.mat'), 'ephysHistoTable')
probes = ephysHistoTable; clear ephysHistoTable
mice = unique(probes.mouseID);

% get probe lines
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
            [probes.mouseID{i} '_registration.mat']), 'tforms');
        load(fullfile(getenv('SSD'), 'paper2', 'histo', 'histoLabels', [probes.mouseID{i} '_histoLabels.mat']), 'scaling', 'ap');
        tform = tforms.(probes.side{i});

        % convert pts to pixels (tform is from histo pixels to ccf pixels)
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
    imgMod = imadjust(img, [0 .8]);

% zoomed out
subplot(1,2,1); hold on
image(d.ml, d.dv, imgMod)
plotLabels2D(d.labels(slice,:,:), 'dvGrid', d.dv, 'mlGrid', d.ml, ...
    'colors', c, 'patchArgs', {'FaceColor', 'none', 'LineWidth', 3, 'LineStyle', '-'});
set(gca, 'Visible', 'off', 'color', 'black')
xlims = xlim;
plot(xlims(2)+[0 -1]-range(xlims)*.05, [5 5], 'Color', 'white', 'linewidth', 3)  % add scale bar
temp = getframe(gca); temp = frame2im(temp); imwrite(temp, 'E:\lab_files\paper2\paper_figures\matlab\histo_zoomout.png')



% zoomed in
subplot(1,2,2); hold on
xlims = [1.2 3.4]; ylims = [1.6 3.6];
image(d.ml, d.dv, imgMod)
plotLabels2D(d.labels(slice,:,:), 'dvGrid', d.dv, 'mlGrid', d.ml, ..., ...
    'colors', c, 'patchArgs', {'FaceColor', 'none', 'LineWidth', 6, 'LineStyle', '-'});
set(gca, 'xlim', xlims, 'ylim', ylims, 'color', 'black')
set(gca, 'Visible', 'off')
plot(xlims(2)+[0 -.5]-range(xlims)*.02, ylims(2)+[0 0]-range(ylims)*.02, 'Color', 'white', 'linewidth', 5)  % plot scale bar
temp = getframe(gca); temp = frame2im(temp); imwrite(temp, 'E:\lab_files\paper2\paper_figures\matlab\histo_zoomin.png')


%% 3D view of all traced stuff

% settings
transparency = .2;
pkcColor = [.8 .6 1];  % purkinje cell layer color



close all
figure('color', 'white', 'position', [2.00 2.00 1278.00 1408.00], 'menubar', 'none'); hold on


bins = strcmp(probes.mouseID, mouse);
x = probes.locHisto(bins,:,1)'; y = probes.locHisto(bins,:,3)';  z = probes.locHisto(bins,:,2)'; % get probe locations

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
    
%     cmat(repmat(squeeze(d.pcLayers(i,:,:) | d.surface(i,:,:)), 1, 1, 3)) = .5;
    cmat(repmat(squeeze(d.surface(i,:,:)), 1, 1, 3)) = 0;
    cmat(repmat(squeeze(d.pcLayers(i,:,:)), 1, 1, 3)) = repelem(pkcColor, sum(d.pcLayers(i,:,:), 'all'));
    
    transp = double(~all(cmat==1, 3)) * transparency/2;
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
                'LineWidth', 1.5, 'EdgeAlpha', .5)
        end
    end
end

% add probe traces
plot3(x, z, y, 'color', [cfg.diiColor], 'LineWidth', 2)  % add probe traces

% scatter units
xyz = cat(1, probes(bins,:).GC_Points{:}) / 1000;
scatter3(xyz(:,1), xyz(:,2), xyz(:,3), 50, [0 0 0], 'filled', 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 1)

% axis props
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

%% show recording locations for all mice!

% this code is absurdly ineffective due to repeated calls to plotLabels2D
% (within plotUnitsOnCcf) // to improve should save 2 dimensional
% projection coordinates from plotLabels2D and plot these manually every
% time to avoid re-computation

% settings
plotProbes = false;
showAllUnits = true;

% figpos = [153.00 226.00 2213.00 1038.00];
figpos = [22.00 178.00 1772.00 1152.00];

close all
figure('color', 'white', 'position', figpos, 'menubar', 'none'); hold on

unitInfo = getUnitInfo('nucleiOnly', ~showAllUnits);

% get colors for scatter points
[~, cind] = ismember(unitInfo.nucleus, {'dentate', 'interpositus', 'fastigial'});
cind(cind==0) = 4;
scatColors = [cfg.nucleusColors; 0 0 0];
scatColors = scatColors(cind,:);

ncols = 4;
nrows = ceil(length(mice) / ncols) * 2;
sind = 1;

for i = 1:length(mice)
    
    bins = strcmp(unitInfo.mouse, mice{i});
    mouseInfo = unitInfo(bins,:);
    
    
    subplots = [nrows ncols sind
                nrows ncols sind+ncols];
    plotUnitsOnCcf(mouseInfo, 'colors', scatColors(bins,:), ...
        'half', false, 'figpos', [2.00 722.00 1278.00 688.00], ...
        'scatArgs', {'MarkerEdgeColor', 'black'}, 'scatSz', 30, ...
        'scalebar', 0, 'subplots', subplots);
    
    if plotProbes; probloc = probes(strcmp(probes.mouseID, mice{i}),:).locCcf; end
    
    subplot(nrows, ncols, sind); hold on
    title(mice{i})
    if plotProbes; lns = plot(probloc(:,:,1)', probloc(:,:,3)', 'LineWidth', 1.5, 'color', cfg.diiColor); uistack(lns, 'bottom'); end
    set(gca, 'visible', 'on', 'xcolor', 'none', 'ycolor', 'none')  % make sure title is visible but not axes
    axpostop = get(gca, 'position');
    
    subplot(nrows, ncols, sind+ncols);hold on
    if plotProbes; lns = plot(probloc(:,:,1)', probloc(:,:,2)', 'LineWidth', 1.5, 'color', cfg.diiColor); uistack(lns, 'bottom'); end
    axposbot = get(gca, 'position');
    vshift = axpostop(2)-axposbot(2)-axposbot(4);  % shift up to butt against plot above
    set(gca, 'position', [axposbot(1) axposbot(2)+vshift axposbot(3) axposbot(4)])
    
    % increment subplot index
    if mod(sind, ncols)==0; sind=sind+ncols; end
    sind=sind+1;
    pause(.1)
end

savefig('E:\lab_files\paper2\paper_figures\matlab\unitlocations_allmice.fig')
saveas(gcf, 'E:\lab_files\paper2\paper_figures\matlab\unitlocations_allmice.svg')

%%
figimg = print('-r600', '-RGBImage');
figure; image(figimg)


















