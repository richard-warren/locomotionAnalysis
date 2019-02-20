function plotOneVsTwoStepTrajectories(data, rowInds, plotType)


% global settings
obsDiam = 3.175; % mm
colors = [.25 1 1; .25 1 .25]; % colors for little vs big steps
obsColor = [.2 .2 .2];
controlColor = [0 0 0];
kinWidthPortion = .8; % portion of fig occupied by kin data
kinematicsHistoSeparation = .02; % separation between kin data dn histos
sidePadding = .02; % padding at side of figure

% kinematic plot settings
xLims = [-.1 .04];
zLims = [0 .015];
circSize = 100;
lineWid = 3;
tracesPerPlot = 5;

% histogram settings
yLimsHisto = [0 .15];
xLimsHisto = [-.06 .06];
xRes = .001;
gausKernelSig = .004; % (m)
transparency = .4;



keyboard

% initializations
% if showPredictedLocations; predictedLocations = [data.swingStartDistance] + [data.predictedLengths]; end % predicted distance to obs
binNum = max(rowInds);
rowInds(~validBins) = 0; % makes sure we don't use invalid trials
deltaLengths = cellfun(@(x,ind) x(1,data(ind).firstModPaw), ... % for each trial, get length of first mod paw
    {data.modifiedSwingLengths}, num2cell(1:length(data))) ...    % all mod paw lengths, trial ind
    - [data.predictedLengths];                                    % subtract predicted lengths
deltaControlLengths = cellfun(@(x,ind) x(end,data(ind).firstModPaw), ...
    {data.controlSwingLengths}, num2cell(1:length(data))) ...
    - [data.predictedControlLengths];
numModSteps = cellfun(@(x,ind) x(1,data(ind).firstModPaw), ...
    {data.modStepNum}, num2cell(1:length(data)));

kernel = arrayfun(@(x) (1/(gausKernelSig*sqrt(2*pi))) * exp(-.5*(x/gausKernelSig)^2), ...
    -gausKernelSig*5:xRes:gausKernelSig*5);
kernel = kernel / sum(kernel);
deltaBinCenters = xLimsHisto(1)-5*gausKernelSig : xRes : xLimsHisto(2)+5*gausKernelSig;
deltaBinEdges = [deltaBinCenters deltaBinCenters(end)+xRes] - .5*xRes;
obsRadius = obsDiam/1000/2;

% figure layout settings
kinematicsWidth = (1-2*sidePadding-kinematicsHistoSeparation) * kinWidthPortion;
histoWidth = (1-2*sidePadding-kinematicsHistoSeparation) * (1-kinWidthPortion);

% close all
figure('color', 'white', 'menubar', 'none', 'position', [100 50 800 1000], 'InvertHardcopy', 'off');

% plot probability density functions (referred to as histos in this script sometimes... argh)
for h = 1:binNum

    ax = subaxis(binNum, 2, h*2);
    oneStepBins = rowInds==h & numModSteps==1;
    twoStepBins = rowInds==h & numModSteps==2;
    oneTwoRatio = sum(oneStepBins) / (sum(oneStepBins) + sum(twoStepBins));

    % control histo
    binCounts = histcounts(deltaControlLengths(rowInds==h), deltaBinEdges);
    histoConv = conv(binCounts, kernel, 'same');
    histoConv = (histoConv/sum(histoConv));
    shadedErrorBar(deltaBinCenters, histoConv, cat(1, histoConv, zeros(1,length(histoConv))), ...
            'lineprops', {'linewidth', 3, 'color', controlColor}, ...
            'patchSaturation', 1-transparency); hold on;
        
    % one step histo
    if any(oneStepBins)
        binCounts = histcounts(deltaLengths(oneStepBins), deltaBinEdges);
        histoConv = conv(binCounts, kernel, 'same');
        histoConv = (histoConv/sum(histoConv)) * oneTwoRatio;
        shadedErrorBar(deltaBinCenters, histoConv, cat(1, histoConv, zeros(1,length(histoConv))), ...
            'lineprops', {'linewidth', 3, 'color', colors(2,:)}, ...
            'patchSaturation', 1-transparency); hold on;
    end

    % two step histo
    binCounts = histcounts(deltaLengths(twoStepBins), deltaBinEdges);
    histoConv = conv(binCounts, kernel, 'same');
    histoConv = (histoConv/sum(histoConv)) * (1-oneTwoRatio);
    shadedErrorBar(deltaBinCenters, histoConv, cat(1, histoConv, zeros(1,length(histoConv))), ...
            'lineprops', {'linewidth', 3, 'color', colors(1,:)}, ...
            'patchSaturation', 1-transparency); hold on;

    % set apearance
    axPos = get(gca, 'position'); axPos(1) = sidePadding+kinematicsWidth+kinematicsHistoSeparation; axPos(3) = histoWidth;
    set(ax, 'box', 'off', 'xlim', xLimsHisto, 'ylim', yLimsHisto, 'tickdir', 'out', ...
        'xtick', [-abs(min(xLimsHisto))*.5 0 abs(min(xLimsHisto))*.5], 'position', axPos, 'ticklength', [0.04 0.025]);
    ax.YAxis.Visible = 'off';
    if h==binNum; xlabel('\Delta swing length (m)'); end
end








% plot trial trajectories
if strcmp(plotType, 'trials')

    for h = 1:binNum

        subaxis(binNum, 2, h*2-1);


        % get inds for one and two step trials within bin
        % (making sure the number of trials of each type reflects the proportion of each type across all trials)
        oneStepPortion = sum(rowInds==h & numModSteps==1) / sum(rowInds==h);
        oneTwoStepTrialNums = [ceil(tracesPerPlot*oneStepPortion) ceil(tracesPerPlot*(1-oneStepPortion))];

        if oneTwoStepTrialNums(1)>0
            oneStepInds = find(rowInds==h & numModSteps==1);
            oneStepInds = oneStepInds(randperm(length(oneStepInds), oneTwoStepTrialNums(1)));
        else
            oneStepInds = [];
        end
        twoStepInds = find(rowInds==h & numModSteps>1);
        twoStepInds = twoStepInds(randperm(length(twoStepInds), oneTwoStepTrialNums(2)));


        % plot individual trials
        for i = [oneStepInds twoStepInds]

            paw = data(i).firstModPaw;
            
            % scatter position of swing foot at obsPos
            scatter(data(i).locations(data(i).obsPosInd,1,paw), ...
                    data(i).locations(data(i).obsPosInd,3,paw), ...
                    50, obsColor, 'filled'); hold on
            
            % plot x and y trajectories
            locationInds = data(i).modifiedStepIdentities(:,paw)==1;
            x = data(i).locations(locationInds,1,paw);
            z = data(i).locations(locationInds,3,paw);

            if ismember(i, oneStepInds); colorInd=2; else; colorInd=1; end
            plot(x, z, 'color', colors(colorInd,:), 'linewidth', 1.5); hold on

            % scatter dots at end of each swing
%             scatter(x(end), z(end), 50, colors(colorInd,:), 'filled'); hold on
        end
        
        % add obstacle
        avgObsHgt = mean([data(rowInds==h).obsHeightsVid]) / 1000; % get avg obstacle height for bin
        circ = rectangle('position', [0-obsRadius, avgObsHgt-2*obsRadius, 2*obsRadius, 2*obsRadius], ...
            'curvature', [1 1], 'facecolor', obsColor, 'edgecolor', 'none');
        
        % add line to bottom
        line(xLims, [0 0], 'color', obsColor, 'linewidth', 2)
        
        % get right control locations
        controlLocations = cellfun(@(x,ind) x{data(ind).firstModPaw}(end,:,:), ...
            {data(rowInds==h).controlLocations}, num2cell(find(rowInds==h)), 'UniformOutput', false); % only take last control step
        controlLocations = cat(1,controlLocations{:});

        % get x offset (mean starting x pos of modified steps)
        modLocations = cellfun(@(x,ind) x{data(ind).firstModPaw}(1,:,:), ...
            {data(rowInds==h).modifiedLocations}, num2cell(find(rowInds==h)), 'UniformOutput', false); % only take last control step
        modLocations = cat(1,modLocations{:});
        xOffset = mean(squeeze(modLocations(:,1,1)));

        x = squeeze(controlLocations(:,1,:));
        x = x - (mean(x(:,1)) - xOffset);
        z = squeeze(controlLocations(:,3,:));
        temp = plot(mean(x,1), mean(z,1), 'color', controlColor, 'linewidth', lineWid); hold on;
        uistack(temp, 'bottom')


        % set appearance
        daspect([1 1 1]);
        axPos = get(gca, 'position'); % axPos(2) = histoHgt+botPadding; axPos(4) = (1-(histoHgt+botPadding)-topPadding);
        axPos(1) = sidePadding; axPos(3) = kinematicsWidth;
        set(gca, 'xlim', xLims, 'ylim', zLims, 'position', axPos);
        axis off
    end
    
    
% plot average trajectories
elseif strcmp(plotType, 'averages')
    
    for h = 1:binNum

        % get subplot bins
        subaxis(binNum, 2, h*2-1);
        
        % get subplot bins for different conditions
        controlBins = rowInds==h;
        oneStepBins = rowInds==h & numModSteps==1;
        twoStepBins = rowInds==h & numModSteps==2;
        oneTwoRatio = sum(oneStepBins) / (sum(oneStepBins) + sum(twoStepBins)); % ratio of trials in which swing foot takes one large step to those in which an additional step is taken
        oneTwoRatio = oneTwoRatio * 2 - 1; % scale from -1 to 1

        
        % plot one step means
        allBins = {twoStepBins, oneStepBins};
        signs = [-1 1]; 

        for i = 1:2
            if any(allBins{i})
                locations = cellfun(@(x,ind) x{data(ind).firstModPaw}(1,:,:), ...
                    {data(allBins{i}).modifiedLocationsInterp}, num2cell(find(allBins{i})), 'UniformOutput', false);
                locations  = cat(1,locations{:});

                x = median(squeeze(locations (:,1,:)),1);
                z = median(squeeze(locations (:,3,:)),1);
                plot(x, z, 'color', colors(i,:), 'linewidth', lineWid + signs(i)*oneTwoRatio*lineWid); hold on;

                contactInds = cellfun(@(x,ind) x(1,data(ind).firstModPaw), ...
                    {data(allBins{i}).pawObsPosIndInterp}, num2cell(find(allBins{i})));
                meanContactInd = round(nanmean(contactInds));
                scatter(x(meanContactInd), z(meanContactInd), ...
                    circSize + signs(i)*oneTwoRatio*circSize, colors(i,:), 'filled'); hold on
            end
        end



        % add obstacle
        avgObsHgt = mean([data(rowInds==h).obsHeightsVid]) / 1000; % get avg obstacle height for bin
        circ = rectangle('position', [0-obsRadius, avgObsHgt-2*obsRadius, 2*obsRadius, 2*obsRadius], ...
            'curvature', [1 1], 'facecolor', obsColor, 'edgecolor', 'none');
        
        % add line to bottom
        line(xLims, [0 0], 'color', obsColor, 'linewidth', 2)
        
        % get control locations
        controlLocations = cellfun(@(x,ind) x{data(ind).firstModPaw}(end,:,:), ...
            {data(rowInds==h).controlLocations}, num2cell(find(rowInds==h)), 'UniformOutput', false); % only take last control step
        controlLocations = cat(1,controlLocations{:});

        % get x offset (mean starting x pos of modified steps)
        modLocations = cellfun(@(x,ind) x{data(ind).firstModPaw}(1,:,:), ...
            {data(rowInds==h).modifiedLocations}, num2cell(find(rowInds==h)), 'UniformOutput', false); % only take last control step
        modLocations = cat(1,modLocations{:});
        xOffset = mean(squeeze(modLocations(:,1,1)));

        x = squeeze(controlLocations(:,1,:));
        x = x - (mean(x(:,1)) - xOffset);
        z = squeeze(controlLocations(:,3,:));
        temp = plot(mean(x,1), mean(z,1), 'color', controlColor, 'linewidth', lineWid); hold on;
        uistack(temp, 'bottom')


        % set appearance
        daspect([1 1 1]);
        axPos = get(gca, 'position');
        axPos(1) = sidePadding; axPos(3) = kinematicsWidth;
        set(gca, 'xlim', xLims, 'ylim', zLims, 'position', axPos);
        axis off
    end
end

blackenFig







