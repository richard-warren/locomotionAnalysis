
unit = 3;

figure;

yyaxis left
plot(obsTimes, obsPositions)

yyaxis right
plot(spkTimes, spkRates(unit,:))

hold on;
scatter(unitTimes{unit}, zeros(1,length(unitTimes{unit})));

pimpFig