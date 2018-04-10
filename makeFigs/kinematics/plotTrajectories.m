function plotTrajectories(data, bins, binLabels, plotType)




% global settings
colors = [.25 1 1; .25 1 .25];
controlColor = repmat(0, 1, 3);
botPadding = .08;
topPadding = .05;
histoHgt = .15;
showPredictedLocations = true;

% kinematic plot settings
yLims = [-.1 .1];
xLims = [-.02 .02];
scaleBarSize = .01;
circSize = 150;
lineWid = 4;
tracesPerPlot = 8;

% histogram settings
yLimsHisto = [0 .1];
xLimsHisto = [-.06 .06];
xRes = .001;
gausKernelSig = .004; % (m)
transparency = .4;



% initializations
if showPredictedLocations
     predictedLocations = [data.swingStartDistance] + [data.predictedLengths]; % predicted distance to obs
end
binNum = max(bins);
deltaLengths = cellfun(@(x) x(1,3), {data.modifiedSwingLengths}) - [data.predictedLengths];
deltaControlLengths = cellfun(@(x) x(2,3), {data.controlSwingLengths}) - [data.predictedControlLengths];
numModSteps = reshape([data.modStepNum],4,length(data))';

kernel = arrayfun(@(x) (1/(gausKernelSig*sqrt(2*pi))) * exp(-.5*(x/gausKernelSig)^2), ...
    -gausKernelSig*5:xRes:gausKernelSig*5);
kernel = kernel / sum(kernel);
deltaBinCenters = min(min(deltaLengths)-2*gausKernelSig, xLimsHisto(1)) : xRes : max(max(deltaLengths)+2*gausKernelSig, xLimsHisto(2));
deltaBinEdges = [deltaBinCenters deltaBinCenters(end)+xRes] - .5*xRes;


close all
figure('color', 'white', 'menubar', 'none', 'position', [100 100 1280 720], 'InvertHardcopy', 'off');

% plot PDFs
for h = 1:binNum

    ax = subaxis(2, binNum , h+binNum, 'marginleft', .01, 'marginright', .01);
    binBins = (bins==h)';
    oneStepBins = binBins & numModSteps(:,3)==1;
    twoStepBins = binBins & numModSteps(:,3)==2;
    oneTwoRatio = sum(oneStepBins) / (sum(oneStepBins) + sum(twoStepBins));

    % control histo
    binCounts = histcounts(deltaControlLengths(binBins), deltaBinEdges);
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
    axPos = get(gca, 'position'); axPos(2) = botPadding; axPos(4) = histoHgt;
    set(ax, 'box', 'off', 'xlim', xLimsHisto, 'ylim', yLimsHisto, 'tickdir', 'out', ...
        'xtick', [-abs(min(xLimsHisto))*.5 0 abs(min(xLimsHisto))*.5], 'position', axPos, 'ticklength', [0.04 0.025]);
    ax.YAxis.Visible = 'off';
    if h==ceil(binNum/2); xlabel('\Delta swing length (m)'); end
end








% plot trial trajectories
if strcmp(plotType, 'trials')

    % initializations
    modStepNum = cellfun(@(x) x(1,3), {data.modStepNum});


    for h = 1:binNum

        subaxis(2, binNum, h, 'marginleft', .01, 'marginright', .01);
        line(xLims, [0 0], 'color', [0 0 0], 'linewidth', 5); hold on; % add line for obstacle
        if showPredictedLocations
            line(xLims, repmat(mean(predictedLocations(bins==h)),1,2), 'color', [0 0 0], 'linewidth', 3, 'linestyle', ':'); hold on; % add line for obstacle
        end


        % get inds for one and two step trials within bin
        % (making sure the number of trials of each type reflects the proportion of each type across all trials)
        binBins = (bins==h);
        oneStepPortion = sum(modStepNum(binBins)==1) / sum(binBins);
        oneTwoStepTrialNums = [round(tracesPerPlot*oneStepPortion) round(tracesPerPlot*(1-oneStepPortion))];

        if oneTwoStepTrialNums(1)>0
            oneStepInds = find(binBins & modStepNum==1);
            oneStepInds = oneStepInds(randperm(length(oneStepInds), oneTwoStepTrialNums(1)));
        else
            oneStepInds = [];
        end
        twoStepInds = find(binBins & modStepNum>1);
        twoStepInds = twoStepInds(randperm(length(twoStepInds), oneTwoStepTrialNums(2)));


        % plot individual trials
        for i = [oneStepInds twoStepInds]
            for j = 2:3

                % plot x and y trajectories
                locationInds = data(i).modifiedStepIdentities(:,j)==1;
                x = data(i).locations(locationInds,1,j);
                y = data(i).locations(locationInds,2,j);
                if data(i).modStepNum(1,j)~=1; colorInd=1; else; colorInd=2; end
                plot(y, x, 'color', colors(colorInd,:), 'linewidth', 1.5); hold on

                % scatter dots at start of each swing
                scatter(y(end), x(end), 100, colors(colorInd,:), 'filled'); hold on

                % scatter position of swing foot at obsPos
                if j==3
                    scatter(data(i).locations(data(i).obsPosInd,2,j), data(i).locations(data(i).obsPosInd,1,j), ...
                        100, colors(colorInd,:), 'x'); hold on
                end
            end
        end


%         % get right control locations
%         controlLocations = cellfun(@(x) x{3}(end,:,:), ...
%             {data(binBins).controlLocationsInterp}, 'uniformoutput', 0); % only take last control step
%         controlLocations = cat(1,controlLocations{:});
% 
%         % get x offset (mean starting x pos of modified steps)
%         modLocations = cellfun(@(x) x{3}(1,:,:), {data(binBins).modifiedLocationsInterp}, 'uniformoutput', 0);
%         modLocations = cat(1,modLocations{:});
%         xOffset = mean(squeeze(modLocations(:,1,1)));
% 
%         x = squeeze(controlLocations(:,1,:));
%         x = x - (mean(x(:,1)) - xOffset);
%         y = squeeze(controlLocations(:,2,:));
%         plot(mean(y,1), mean(x,1), 'color', controlColor, 'linewidth', lineWid); hold on;




        % set appearance
        axPos = get(gca, 'position'); axPos(2) = histoHgt+botPadding; axPos(4) = (1-(histoHgt+botPadding)-topPadding);
        set(gca, 'dataaspectratio', [1 1 1], 'xlim', [-.02 .02], 'ylim', yLims, 'position', axPos);
        axis off
    %     xlabel(['predicted dist to obs (m): ' binLabels{h}]);
    end
    
    
% plot average trajectories
elseif strcmp(plotType, 'averages')
    
    % initializations
    obsPosIndInterps = cellfun(@(x) x(1,3), {data.pawObsPosIndInterp});


    for h = 1:binNum

        % get subplot bins
        subaxis(2, binNum, h, 'marginleft', .01, 'marginright', .01);
        line(xLims, [0 0], 'color', [0 0 0], 'linewidth', 5); hold on; % add line for obstacle
        if showPredictedLocations
            line(xLims, repmat(mean(predictedLocations(bins==h)),1,2), 'color', [0 0 0], 'linewidth', 3, 'linestyle', ':'); hold on; % add line for obstacle
        end
        
        % get subplot bins for different conditions
        binBins = (bins==h)';
        controlBins = binBins;
        leftModBins = binBins & numModSteps(:,2)==1;
        rightModOneStepBins = binBins & (numModSteps(:,3)==1);
        rightModTwoStepBins = binBins & (numModSteps(:,3)==2);
        oneTwoRatio = sum(rightModOneStepBins) / (sum(rightModOneStepBins) + sum(rightModTwoStepBins)); % ratio of trials in which swing foot takes one large step to those in which an additional step is taken
        oneTwoRatio = oneTwoRatio * 2 - 1; % scale from -1 to 1


        % get left and right control locations
        controlLocations = {data(controlBins).controlLocationsInterp};
        leftControlLocations = cellfun(@(x) x{2}(end,:,:), controlLocations, 'uniformoutput', 0);
        leftControlLocations = cat(1,leftControlLocations{:});
        rightControlLocations = cellfun(@(x) x{3}(end,:,:), controlLocations, 'uniformoutput', 0);
        rightControlLocations = cat(1,rightControlLocations{:});

        % get left modified locations
        leftModLocations = {data(leftModBins).modifiedLocationsInterp};
        leftModLocations = cellfun(@(x) x{2}, leftModLocations, 'uniformoutput', 0);
        leftModLocations = cat(1,leftModLocations{:});

        % get right modified (one step) locations
        rightModOneStepLocations = {data(rightModOneStepBins).modifiedLocationsInterp};
        rightModOneStepLocations = cellfun(@(x) x{3}, rightModOneStepLocations, 'uniformoutput', 0);
        rightModOneStepLocations = cat(1,rightModOneStepLocations{:});

        % get right modified (two step) locations
        rightModTwoStepLocations = {data(rightModTwoStepBins).modifiedLocationsInterp};
        rightModTwoStepLocations = cellfun(@(x) x{3}(1,:,:), rightModTwoStepLocations, 'uniformoutput', 0);
        rightModTwoStepLocations = cat(1,rightModTwoStepLocations{:});

        % get left and right x offsets
        modLocations = cellfun(@(x) x{2}(1,:,:), {data(binBins).modifiedLocationsInterp}, 'uniformoutput', 0);
        modLocations = cat(1,modLocations{:});
        xOffsetLeft = mean(squeeze(modLocations(:,1,1)));
        modLocations = cellfun(@(x) x{3}(1,:,:), {data(binBins).modifiedLocationsInterp}, 'uniformoutput', 0);
        modLocations = cat(1,modLocations{:});
        xOffsetRight = mean(squeeze(modLocations(:,1,1)));


%         % plot control left
%         x = squeeze(leftControlLocations(:,1,:));
%         x = x - (mean(x(:,1)) - xOffsetLeft);
%         y = squeeze(leftControlLocations(:,2,:));
%         plot(mean(y,1), mean(x,1), 'color', controlColor, 'linewidth', lineWid); hold on;
% 
%         % plot control right
%         x = squeeze(rightControlLocations(:,1,:));
%         x = x - (mean(x(:,1)) - xOffsetRight);
%         y = squeeze(rightControlLocations(:,2,:));
%         plot(mean(y,1), mean(x,1), 'color', controlColor, 'linewidth', lineWid); hold on;

        % plot mod left
        x = squeeze(leftModLocations(:,1,:));
        y = squeeze(leftModLocations(:,2,:));
        plot(mean(y,1), mean(x,1), 'color', colors(2,:), 'linewidth', lineWid); hold on;

        % plot mod right, one step
        if ~isempty(rightModOneStepLocations)
            x = squeeze(rightModOneStepLocations(:,1,:));
            y = squeeze(rightModOneStepLocations(:,2,:));
            plot(mean(y,1), mean(x,1), 'color', colors(2,:), 'linewidth', lineWid + oneTwoRatio*lineWid); hold on;
        end

        % plot mod right, two step
        if ~isempty(rightModTwoStepLocations)
            x = squeeze(rightModTwoStepLocations(:,1,:));
            y = squeeze(rightModTwoStepLocations(:,2,:));
            plot(mean(y,1), mean(x,1), 'color', colors(1,:), 'linewidth', lineWid + -oneTwoRatio*lineWid); hold on;
        end

        % mark avg position of obsPos
        if any(rightModOneStepBins)
            oneStepObsPos = round(mean(obsPosIndInterps(rightModOneStepBins)));
            scatter(mean(squeeze(rightModOneStepLocations(:,2,oneStepObsPos))), ...
                    mean(squeeze(rightModOneStepLocations(:,1,oneStepObsPos))), circSize + oneTwoRatio*circSize, colors(2,:), 'filled'); hold on
        end

        if any(rightModTwoStepBins)
            twoStepObsPos = round(nanmean(obsPosIndInterps(rightModTwoStepBins)));
            scatter(mean(squeeze(rightModTwoStepLocations(:,2,twoStepObsPos))), ...
                    mean(squeeze(rightModTwoStepLocations(:,1,twoStepObsPos))), circSize + -oneTwoRatio*circSize, colors(1,:), 'filled'); hold on
        end



        % set appearance
        axPos = get(gca, 'position'); axPos(2) = histoHgt+botPadding; axPos(4) = (1-(histoHgt+botPadding)-topPadding);
        set(gca, 'dataaspectratio', [1 1 1], 'xlim', [-.02 .02], 'ylim', yLims, 'position', axPos);
        axis off
%         xlabel(['predicted dist to obs (m): ' binLabels{h}]);
    end
end






% add scale bar
line([xLims(2)-scaleBarSize xLims(2)], repmat(yLims(2),1,2)-.02, 'linewidth', 3, 'color', 'black')
text(xLims(2)-.5*scaleBarSize, yLims(2)-.02+.005, sprintf('%i mm', scaleBarSize*1000), 'horizontalalignment', 'center')

saveas(gcf, [getenv('OBSDATADIR') 'figures\trialKinematics.png']);
savefig([getenv('OBSDATADIR') 'figures\trialKinematics.fig'])
blackenFig;
print('-clipboard', '-dmeta')






