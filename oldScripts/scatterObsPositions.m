


swingLengths = reshape([data.modifiedSwingLengths]',4,[])'; swingLengths = swingLengths(:,3)';
vels = reshape([data.modifiedWheelVels]',4,[])'; vels = vels(:,3)';
predictedDistToObs = [data.swingStartDistance] + [data.predictedLengths];
actualDistToObs = [data.swingStartDistance] + swingLengths;

close all; figure('menubar', 'none', 'color', 'white', 'position', [645   329   777   639]);
colors = winter(length(data));
[~,sortInds] = sort(vels);

scatter(predictedDistToObs(sortInds)*1000, actualDistToObs(sortInds)*1000, 50, colors, 'filled');
daspect([1 1 1])
xlabel('predicted distance to obs (mm)')
ylabel('actual distance to obs (mm)')

saveas(gcf, [getenv('OBSDATADIR') 'figures\distToObsScatter.png']);


figure('menubar', 'none', 'color', 'white', 'position', [645   329   777   639]);
numModSteps = reshape([data.modStepNum],4,length(data))';
shortStepBins = numModSteps(:,3)>1;

histogram(predictedDistToObs(shortStepBins)*1000, 'binwidth', 2.5, 'facecolor', [221 0 239]/255); hold on;
histogram(actualDistToObs(shortStepBins)*1000, 'binwidth', 2.5, 'facecolor', [0 0 0]);
set(gca, 'ycolor', 'white', 'box', 'off')

saveas(gcf, [getenv('OBSDATADIR') 'figures\distToObsHisto.png']);