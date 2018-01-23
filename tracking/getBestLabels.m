function bestLabels = getBestLabels(scores, objectNum, isOccluded)

% !!! need to document


sequences = perms(1:objectNum);
labels = nan(size(sequences));
labelScores = nan(size(sequences));


% first see whether the labels are unique for each object
% if so, simply return these labels and skip the algorithm below
[maxes, bestLabels] = max(scores, [], 2);
if length(unique(bestLabels)) == length(bestLabels)
    bestLabels(maxes==0) = nan; 
    return;
end


% only test sequences in which occluded states are LAST
if any(isOccluded)
    occludedObs = find(isOccluded);
    finalObsInSeq = sequences(:, (end-length(occludedObs))+1:end); % these are the columns corresponding to the final objects in the sequence (ie last column if one occluded object, last two columns if two occluded objects)
    validSequences = all( ismember(finalObsInSeq, repmat(occludedObs, size(sequences,1), 1)), 2); % these are the sequences in which the occluded states are LAST
    sequences = sequences(validSequences, :);
end


% iterate through all possible sequences of objects
for i = 1:size(sequences,1)
    
    scoresTemp = scores;
    
    for j = 1:objectNum
        
        % get index of current object in sequence
        obInd = sequences(i,j);

        % get max score and label for max score
        [maxTemp, labelTemp] = max(scoresTemp(obInd,:));
        labelScores(i, obInd) = maxTemp;
        
        % if max score is zero set label to nan
        if maxTemp>0
            labels(i, obInd) = labelTemp;
            scoresTemp(:, labelTemp) = 0; % set scores for this label to zero for all objects to avoid reassigning this label to another object
        else
            labels(i, obInd) = nan;
        end
    end
end

[~, bestSequence] = max(sum(labelScores,2));
bestLabels = labels(bestSequence, :)';


