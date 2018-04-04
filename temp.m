


load([getenv('OBSDATADIR') 'kinematicData.mat'], 'data')





%% test arduino constant acceleration ramp

a = 5.0;
maxV = 1.5;
% getV2 = @(v1) .5*(v1 + sqrt(4*a+v1^2));
getV2 = @(v1) .5 * (v1 + sqrt(4*a*mPerTick + v1^2));


% initializations
microStepping = 16;
motorSteps = 200;
timingPulleyRad = 15.2789;
mPerTick = ((2*pi*timingPulleyRad) / (1000*(motorSteps*microStepping)));
v2 = 0;
v = [0];
t = [0];

while v<maxV
    
    v1 = v2;
    v2 = getV2(v1);
    v(end+1) = v2;
    t(end+1) = t(end) + (1/v2)/s;
    disp(v2)
    
end

% close all; scatter(t, v); pimpFig
length(v)

%%
controls = {data.controlLocations};
controls = cellfun(@(x) x{3}, controls, 'uniformoutput', 0);
controls = cat(1,controls{:});

close all; figure; plot(mean(squeeze(controls(:,1,:)),1), mean(squeeze(controls(:,2,:)),1))
set(gca, 'dataaspectratio', [1 1 1])


%%

isInSwing = nan(length(data),2);

for i = 1:length(data)
    
    ind = data(i).obsPosInd;
    isLeftSwing = ~isnan(data(i).modifiedStepIdentities(ind,2));
    isRightSwing = ~isnan(data(i).modifiedStepIdentities(ind,3));
    
    isInSwing(i,:) = [isLeftSwing isRightSwing];
    
end

bothSwingInds = find(and(isInSwing(:,1), isInSwing(:,2)));
bothStanceInds = find(sum(isInSwing,2)==0);

%%
oneSwingOneStance = find([data.oneSwingOneStance]);
trial = oneSwingOneStance(3);
trial=6;


ind = data(i).obsPosInd;
close all; figure;
colors = winter(2);
forepaws = 2:3;

for i = 1:length(forepaws)
    
    plot(data(trial).timeStamps, ...
        data(trial).locations(:,1,forepaws(i)), 'linewidth', 2, 'color', colors(i,:)); hold on;
    
    realBins = ~isnan(data(trial).controlStepIdentities(:,forepaws(i)));
    for j = unique(data(trial).controlStepIdentities(realBins,forepaws(i)))'
        bins = data(trial).controlStepIdentities(:,forepaws(i)) == j;
        plot(data(trial).timeStamps(bins), ...
            data(trial).locations(bins,1,forepaws(i)), 'linewidth', 10, 'color', [0 0 0]); hold on;
    end
    
    realBins = ~isnan(data(trial).modifiedStepIdentities(:,forepaws(i)));
    for j = unique(data(trial).modifiedStepIdentities(realBins,forepaws(i)))'
        bins = data(trial).modifiedStepIdentities(:,forepaws(i)) == j;
        plot(data(trial).timeStamps(bins), ...
            data(trial).locations(bins,1,forepaws(i)), 'linewidth', 10, 'color', colors(i,:)); hold on;
    end
    
end

% line([ind ind], get(gca,'ylim'))
xlims = get(gca,'xlim');
line(xlims, [0 0])
line(repmat(data(trial).timeStamps(data(trial).obsPosInd),1,2), get(gca,'ylim'))
pimpFig; set(gca,'xlim',xlims)


%%

close all; figure;

plot(data(i).locations(:,2,2)); hold on; plot(data(i).locations(:,2,3))







