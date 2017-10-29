% load tracking data
load('C:\Users\LindseyBuckingham\Google Drive\columbia\obstacleData\svm\trackedData\tracked.mat', 'locations')

% user settings
maxVelocity = 15;
velocityWeight = 1;

% initializations
frameHeight = 242; % !!! hacky temp
frameWidth = 398;



% compute unary potentials
unaries = cell(1,length(locations));

for i = 1:length(locations)
    
    frameUnaries = nan(4, length(locations(i).x));
    
    frameUnaries(1,:) = getUnaryPotentials(locations(i).x, locations(i).y, frameWidth, frameHeight, 0, 0); % LH
    frameUnaries(2,:) = getUnaryPotentials(locations(i).x, locations(i).y, frameWidth, frameHeight, 0, 1); % LF
    frameUnaries(3,:) = getUnaryPotentials(locations(i).x, locations(i).y, frameWidth, frameHeight, 1, 0); % RH
    frameUnaries(4,:) = getUnaryPotentials(locations(i).x, locations(i).y, frameWidth, frameHeight, 1, 1); % RF
    
    unaries{i} = frameUnaries;
    
end



%% compute pairwise potentials
pairwise = cell(1,length(locations)-1);

for i = 1:length(locations)-1
    
    framePairwise = getPairwisePotentials([locations(i).x, locations(i).y], [locations(i+1).x, locations(i+1).y], maxVelocity, velocityWeight);
    pairwise{i} = framePairwise;
    
end
