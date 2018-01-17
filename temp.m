

load('C:\Users\rick\Google Drive\columbia\obstacleData\sessions\171202_000\tracking\runTopHandLabeledLocationsObsFrames.mat')
locationFrameIndsTemp = locationFrameInds;
locationsTemp = locations;

load('C:\Users\rick\Google Drive\columbia\obstacleData\sessions\171202_000\tracking\runTopHandLabeledLocationsObsOnFrames.mat')
locationFrameIndsTemp = [locationFrameIndsTemp locationFrameInds];
locationsTemp = cat(2, locationsTemp, locations);

nanBins = isnan(locationFrameIndsTemp);
locations = locationsTemp(~nanBins);
locationFrameInds = locationFrameIndsTemp(~nanBins);

clear locationFrameIndsTemp locationsTemp nanBins
save('C:\Users\rick\Google Drive\columbia\obstacleData\sessions\171202_000\tracking\runTopHandLabeledLocations.mat',...
    'locations', 'locationFrameInds')

%%

I = eye(max(labels));
labelsHotOne = I(labels,:)';