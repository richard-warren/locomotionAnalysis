function splitTrials = splitByRewards(data, dataTimes, rewardTimes, normalize)
    
    % splits continuosly recorded data (e.g. wheel position, velocity) into a cell array where each entry contains only the data recorded before a given reward
    % the number of entires will equal the number of rewards
    % important: this code skips the first trial
    %
    % input        data:          continuosly recorded data
    %              dataTimes:     timestamps for data
    %              rewardTimes:   times of rewards
    %              normalize:     boolean // if 1, the starting value of each trial is subtracted from each trial (used to normalize positional values - don't use for velocity!)
    %
    % output       splitTrials:   cell array where  each entry contains only the data recorded before a given reward


    % initializations
    splitTrials = cell(size(rewardTimes));

    
    % iterate through all rewards
    for i = 2:length(rewardTimes)

        % extract data for single trial
        trialData = data(dataTimes>rewardTimes(i-1) & dataTimes<rewardTimes(i));
    
        % normalize
        if normalize
            trialData = trialData - trialData(1);
        end

        % store data
        splitTrials{i-1} = trialData;
    end
end





% end