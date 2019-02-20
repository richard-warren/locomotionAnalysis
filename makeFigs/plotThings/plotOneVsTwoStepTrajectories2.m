function plotOneVsTwoStepTrajectories2(data, rowInds)

% to do: add individual trials option // add times of whisker contact


% temp
% data = getNestedStructFields(data, {'mouse', 'session', 'trial', 'modPawKinInterp', 'preModPawKinInterp', 'isBigStep', ...
%     'modPawPredictedDistanceToObs', 'preModPawDeltaLength', 'modPawDeltaLength', 'obsHgt'});
% lims = prctile([data.modPawPredictedDistanceToObs], [5 95]);
% rowInds = discretize([data.modPawPredictedDistanceToObs], linspace(lims(1), lims(2), 4));


% global settings
obsRadius = 3.175/1000/2; % (m)
colors = [.25 1 1; .25 1 .25]; % colors for little vs big steps
obsColor = [.2 .2 .2];
controlColor = [0 0 0];
kinWidthPortion = .8; % portion of fig occupied by kin data
kinematicsHistoSeparation = .02; % separation between kin data and histos
sidePadding = .02; % padding at side of figure

% kinematic plot settings
xLims = [-.08 .03];
zLims = [0 .015];
lineWid = 3;

% histogram settings
histLims = [-.04 .04];
bins = 100;
transparency = .4;



% initializations
numRows = max(rowInds);
kinematicsWidth = (1-2*sidePadding-kinematicsHistoSeparation) * kinWidthPortion;
histoWidth = (1-2*sidePadding-kinematicsHistoSeparation) * (1-kinWidthPortion);
plotHgt = (diff(zLims)/diff(xLims))*kinematicsWidth;
histGrid = linspace(histLims(1), histLims(2), bins);
figure('color', 'white', 'menubar', 'none', 'position', [2000 0 800 200*numRows], 'InvertHardcopy', 'off');




% plot probability density functions (referred to as histos in this script sometimes... argh)
for h = 1:numRows

    ax = subplot(numRows, 2, h*2);
    oneStepBins = rowInds==h & [data.isBigStep];
    twoStepBins = rowInds==h & ~[data.isBigStep];
    oneTwoRatio = sum(oneStepBins) / (sum(oneStepBins) + sum(twoStepBins));

    % control histo
    histoConv = ksdensity([data(rowInds==h).preModPawDeltaLength], histGrid);
    shadedErrorBar(histGrid, histoConv, cat(1, histoConv, zeros(1,length(histoConv))), ...
            'lineprops', {'linewidth', 3, 'color', controlColor}, ...
            'patchSaturation', 1-transparency); hold on;
    
    % one step histo
    if any(oneStepBins)
        histoConv = ksdensity([data(oneStepBins).modPawDeltaLength], histGrid) * oneTwoRatio;
        shadedErrorBar(histGrid, histoConv, cat(1, histoConv, zeros(1,length(histoConv))), ...
            'lineprops', {'linewidth', 3, 'color', colors(2,:)}, ...
            'patchSaturation', 1-transparency); hold on;
    end

    % two step histo
    histoConv = ksdensity([data(twoStepBins).modPawDeltaLength], histGrid) * (1-oneTwoRatio);
    shadedErrorBar(histGrid, histoConv, cat(1, histoConv, zeros(1,length(histoConv))), ...
            'lineprops', {'linewidth', 3, 'color', colors(1,:)}, ...
            'patchSaturation', 1-transparency); hold on;

    % set apearance
    axPos = [sidePadding+kinematicsWidth+kinematicsHistoSeparation, ...
             1 - (h/(numRows+1)), ...
             histoWidth, ...
             plotHgt];
    set(ax, 'XLim', histLims, 'TickDir', 'out', 'Position', axPos, 'YColor', 'none')
    if h==numRows; xlabel('\Delta swing length (m)'); end
end






for h = 1:numRows

    % get subplot bins
    subplot(numRows, 2, h*2-1);

    % get subplot bins for different conditions
    controlBins = rowInds==h;
    oneStepBins = rowInds==h & [data.isBigStep];
    twoStepBins = rowInds==h & ~[data.isBigStep];
    oneTwoRatio = sum(oneStepBins) / (sum(oneStepBins) + sum(twoStepBins)); % ratio of trials in which swing foot takes one large step to those in which an additional step is taken
    oneTwoRatio = oneTwoRatio * 2 - 1; % scale from -1 to 1


    % plot one step means
    allBins = {twoStepBins, oneStepBins};
    signs = [-1 1]; 

    for i = 1:2
        if any(allBins{i})

            locations = permute(cat(3, data(allBins{i}).modPawKinInterp), [3,1,2]);
            x = nanmedian(squeeze(locations(:,1,:)), 1);
            z = nanmedian(squeeze(locations(:,3,:)), 1);
            plot(x, z, 'color', colors(i,:), 'linewidth', lineWid + signs(i)*oneTwoRatio*lineWid); hold on;

%                 contactInds = cellfun(@(x,ind) x(1,flat(ind).firstModPaw), ...
%                     {flat(allBins{i}).pawObsPosIndInterp}, num2cell(find(allBins{i})));
%                 meanContactInd = round(nanmean(contactInds));
%                 scatter(x(meanContactInd), z(meanContactInd), ...
%                     circSize + signs(i)*oneTwoRatio*circSize, colors(i,:), 'filled'); hold on
        end
    end

    % add obstacle
    avgObsHgt = nanmean([data(rowInds==h).obsHgt]); % get avg obstacle height for bin
    rectangle('position', [0-obsRadius, avgObsHgt-2*obsRadius, 2*obsRadius, 2*obsRadius], ...
        'curvature', [1 1], 'facecolor', obsColor, 'edgecolor', 'none');

    % add line to bottom
    line(xLims, [0 0], 'color', obsColor, 'linewidth', 2)

    % get control locations
    controlLocations = permute(cat(3, data(controlBins).preModPawKinInterp), [3,1,2]);

    % get x offset (mean starting x pos of modified steps)
    modLocations = permute(cat(3, data(oneStepBins | twoStepBins).modPawKinInterp), [3,1,2]);
    xOffset = nanmean(squeeze(modLocations(:,1,1)));

    x = squeeze(controlLocations(:,1,:));
    x = x - (nanmean(x(:,1)) - xOffset);
    z = squeeze(controlLocations(:,3,:));
    temp = plot(nanmedian(x,1), nanmedian(z,1), 'color', controlColor, 'linewidth', lineWid); hold on;
    uistack(temp, 'bottom')


    % set appearance
    set(gca, 'xlim', xLims, 'ylim', zLims, 'DataAspectRatio', [1 1 1]);
    axPos = [sidePadding, ...
             1 - (h/(numRows+1)), ...
             kinematicsWidth, ...
             plotHgt];
    set(gca, 'position', axPos, 'XLim', xLims, 'YLim', zLims, 'visible', 'off')
end






