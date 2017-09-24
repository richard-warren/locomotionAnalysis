file = 'C:\Users\Rick\Google Drive\columbia\analysis\run.mat';

load(file);

% decoder rotary encoder

% get raw data
aTimes = encoderA.times';
aStates = logical(encoderA.level)';
bTimes = encoderB.times';
bStates = logical(encoderB.level)';

% pad beginnings with initial pin states
% aStates = [~aStates(1) aStates];
% aTimes = [0 aTimes];
% bStates = [~bStates(1) bStates];
% bTimes = [0 bTimes];
%

tics = nan([2, length(aTimes) + length(bTimes) + 1]); % add one for the initial state
tics(:,1) = [~aStates(1); ~bStates(1)];               % initialize baseline state of encoders

% sort all time indices
allTimes = [aTimes bTimes];
[allTimes, timesSortedInds] = sort(allTimes);

% create 2xN matrix of encoder states, sorted with time indices determind above
% each row is an encoder (A or B), and NaN values temporarily stand in for times when the other encoder changes state
allStates = nan(2, length(aTimes)+length(bTimes));

temp = [aStates, nan(1, length(bStates))];
allStates(1,:) = temp(timesSortedInds);

temp = [nan(1, length(aStates)), bStates];
allStates(2,:) = temp(timesSortedInds);

%% replace all NaN values with OPPOSITE of subsequent real number

% things to do:
% -convert to real world units
% -get the initial state of the pins so you're not missing the first index
% - combine both for loops into one?

tic
for i = 1:2
    for j = find(isnan(allStates(i,:)))
        
        % find next non-NaN value and store the opposite of that value in current index
        nextValInd = find(~isnan(allStates(i,j:end)), 1, 'first');
        nextVal = ~allStates(i, j+nextValInd-1);
        
        if ~isempty(nextVal)
            allStates(i,j) = nextVal;
        else
            % fill in missing NaN values at the very end
            allStates(i,j:end) = allStates(i,j-1);
        end
    end    
end



% convert state transitions into code for forward (1) or backward (-1) tics
transitions = nan(1, length(allTimes)-1);
positions = nan(1, length(allTimes)-1);
lookUp = [0,-1,1,0,1,0,0,-1,-1,0,0,1,0,1,-1,0] * -1;

currentPos = 0;

for i = 2:length(allStates)
    
    transitionCode = [num2str(allStates(1,i-1))...
                      num2str(allStates(2,i-1))...
                      num2str(allStates(1,i))...
                      num2str(allStates(2,i))];
    
    transitions(i-1) = lookUp(bin2dec(transitionCode));
    currentPos =  currentPos + transitions(i-1);
    positions(i-1) = currentPos;
    
    
end
toc
%%

close all; figure; plot(allTimes(2:end), positions); pimpFig












