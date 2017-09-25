% things to do:
% -combine both for loops into one?


file = 'C:\Users\Rick\Google Drive\columbia\analysis\run.mat';

% wheel, encoder settings
encoderSteps = 2880; % 720cpr * 4
wheelRad = 95.25; % mm
mmPerTic = (2*wheelRad*pi) / encoderSteps;

load(file);

% decoder rotary encoder

% get raw data
aTimes = encoderA.times';
aStates = logical(encoderA.level)';
bTimes = encoderB.times';
bStates = logical(encoderB.level)';

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
positions = nan(1, length(allTimes));
lookUp = [0,-1,1,0,1,0,0,-1,-1,0,0,1,0,1,-1,0] * -1;

% compute first transition
transitionCode = [num2str(~aStates(1))... % the initial state is always the opposite of the state detected at the first event
                  num2str(~bStates(1))...
                  num2str(allStates(1,1))...
                  num2str(allStates(2,1))];
positions(1) = lookUp(bin2dec(transitionCode)); 

% compute remaining transitions
for i = 2:length(allStates)
    
    transitionCode = [num2str(allStates(1,i-1))...
                      num2str(allStates(2,i-1))...
                      num2str(allStates(1,i))...
                      num2str(allStates(2,i))];
    
    positions(i) = positions(i-1) + lookUp(bin2dec(transitionCode)); 
end
toc

% convert to real-world units
positions = positions * mmPerTic / 1000;

%%

close all; figure; plot(allTimes, positions); pimpFig












