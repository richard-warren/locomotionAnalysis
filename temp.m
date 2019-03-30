

binEdges = -.5:20.5;

flat = getNestedStructFields(data, {'mouse', 'session', 'frameCounts', 'touchFrames', 'condition'});

flat = flat(ismember({flat.mouse}, {'sen7', 'sen8'}) & ...
            ismember({flat.condition}, {'pre', 'noWisk'}));

        
%%
close all; figure
histogram([flat(strcmp({flat.condition}, 'pre')).touchFrames], binEdges, 'Normalization', 'probability'); hold on
histogram([flat(strcmp({flat.condition}, 'noWisk')).touchFrames], binEdges, 'Normalization', 'probability');
legend({'pre', 'noWisk'})

%%
thresholds = 0:20;

successes = nan(length(thresholds),2);
for i = 1:length(thresholds)
    successes(i,1) = mean([flat(strcmp({flat.condition}, 'pre')).touchFrames]<thresholds(i));
    successes(i,2) = mean([flat(strcmp({flat.condition}, 'noWisk')).touchFrames]<thresholds(i));
end