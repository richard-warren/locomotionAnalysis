load('C:\Users\rick\Google Drive\columbia\obstacleData\sessions\171202_000\tracking\runBotHandLabeledLocations200 (all obs frames).mat')
locationFrameInds2 = locationFrameInds;
locations2 = locations;
load('C:\Users\rick\Google Drive\columbia\obstacleData\sessions\171202_000\tracking\runBotHandLabeledLocations325 (few obs frames).mat')

locations = cat(2, locations, locations2);
locationFrameInds = [locationFrameInds, locationFrameInds2];

validInds = ~isnan(locationFrameInds);
locations = locations(:, validInds, :);
locationFrameInds = locationFrameInds(validInds);

clear locations2 locationFrameInds2 validInds

save('C:\Users\rick\Google Drive\columbia\obstacleData\sessions\171202_000\tracking\runBotHandLabeledLocations.mat',...
     'locations', 'locationFrameInds')
 
%%

