% load tracking data
load('C:\Users\LindseyBuckingham\Google Drive\columbia\obstacleData\svm\trackedData\tracked.mat', 'locations')

% user settings
occludedCost = .01;
occludedStates = 1;
maxVelocity = 15;
velocityWeight = 1;

% initializations
frameHeight = 242; % !!! hacky temp
frameWidth = 398;



% compute unary potentials
unary = cell(1,length(locations));

for i = 1:length(locations)
    
    frameUnaries = nan(4, length(locations(i).x) + occludedStates);
    
    frameUnaries(1,1:end-occludedStates) = getUnaryPotentials(locations(i).x, locations(i).y, frameWidth, frameHeight, 0, 0); % LH
    frameUnaries(2,1:end-occludedStates) = getUnaryPotentials(locations(i).x, locations(i).y, frameWidth, frameHeight, 0, 1); % LF
    frameUnaries(3,1:end-occludedStates) = getUnaryPotentials(locations(i).x, locations(i).y, frameWidth, frameHeight, 1, 0); % RH
    frameUnaries(4,1:end-occludedStates) = getUnaryPotentials(locations(i).x, locations(i).y, frameWidth, frameHeight, 1, 1); % RF
    
    frameUnaries(:,end-occludedStates+1:end) = 0;%occludedCost;
    
    unary{i} = frameUnaries; % documentation on match2nd appears to be incorrect... this matrix should be flipped i think...
    
end



% compute pairwise potentials
pairwise = cell(1,length(locations)-1);

for i = 1:length(locations)-1
    
    pairwise{i} = getPairwisePotentials([locations(i+1).x, locations(i+1).y], [locations(i).x, locations(i).y], maxVelocity, velocityWeight, occludedCost);
    
end



%% perform tracking!!! wtf???!!!




% [labels] = match2nd(unary, pairwise, [], 1, 0);
% [labels] = match2nd(cellfun(@(x) x', unaries, 'un', 0), pairwise, [], 1, 0);













