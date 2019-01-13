function [matchedBins, weights] = findMatchedBins(data, conditions, velTolerance, angleTolerance)

% given kinematic or speedAvoidance data struct, finds bins for
% subdistribution within the data that are matched for 'conditions' in
% terms of velocity and body angle // velTolerance and angleTolerance is
% how close in vel and angle two trials have to be to be considered
% 'matched'
% TO DO: perhaps make this more general, s.t. arbitrary variables of
% interest can be matched, and arbitrary number of variables...


% get vel and angle bins
velBinEdges = 0:velTolerance:max([data.avgVel]);
[~, velBinEdges, velBins] = histcounts([data.avgVel], velBinEdges);
angleBinEdges = min(([data.avgAngle])):angleTolerance:max(([data.avgAngle]));
[~, angleBinEdges, angleBins] = histcounts(([data.avgAngle]), angleBinEdges);

controlBins = strcmp({data.condition}, conditions{1});
manipBins = strcmp({data.condition}, conditions{2});
matchedBins = false(1,length(data));
weights = zeros(1,length(data));
mice = unique({data.mouse});


for i = 1:length(velBinEdges)-1
    for j = 1:length(angleBinEdges)-1
        for k = 1:length(mice)

            mouseBins = strcmp({data.mouse}, mice{k});
            binControlTrials = find(velBins==i & angleBins==j & controlBins & mouseBins);
            binManipTrials   = find(velBins==i & angleBins==j & manipBins & mouseBins);

            if ~isempty(binControlTrials) && ~isempty(binManipTrials)
                minCount = min(length(binControlTrials), length(binManipTrials));
                
                matchedBins(binControlTrials(randperm(length(binControlTrials), minCount))) = true;
                matchedBins(binManipTrials(randperm(length(binManipTrials), minCount))) = true;
                
                weights(binControlTrials) = minCount / length(binControlTrials);
                weights(binManipTrials) = minCount / length(binManipTrials);
            end
        end
    end
end

fprintf('%.2f trials in matched sub-distribution\n', mean(matchedBins))