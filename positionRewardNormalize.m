function posNorm = positionRewardNormalize(positions, positionTimes, rewardTimes)

    % normalizes wheel positions such that the positions reset after every reward is reached
    
    posNorm = positions;
    
    for i = 1:length(rewardTimes)
        
        inds = positionTimes>rewardTimes(i);
        rewardPos = posNorm(find(positionTimes<rewardTimes(i), 1, 'last'));
        posNorm(inds) = posNorm(inds) - rewardPos;
        
    end
end