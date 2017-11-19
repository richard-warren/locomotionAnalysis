function getXAlignment


xData = dlmread('xAlignment\bonsai\xTopBot.csv');
xTop = xData(:,1);
xBot = xData(:,2);
validInds = ~isnan(xTop) & ~isnan(xBot);

xLinearMapping = polyfit(xBot(validInds), xTop(validInds), 1);

save('xAlignment\xLinearMapping.mat', 'xLinearMapping')