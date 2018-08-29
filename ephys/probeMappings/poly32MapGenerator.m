neuro = [23 25 27 29 31 19 17 21 11 15 13 1 3 5 7 9;
         24 26 28 30 32 20 18 22 12 16 14 2 4 6 8 10];

intan = [23 22 21 20 19 18 17 16 15 14 13 12 11 10 9 8;
         24 25 26 27 28 29 30 31 0 1 2 3 4 5 6 7];
intan = intan+1;
intan = fliplr(intan);


connectorInds = [16 6 5 15 4 7 3 8 2 9 1 10 14 13 12 11 22 21 20 19 23 25 24 18 26 17 27 29 28 31 30 32];
kiloSortInds = [31 28 25 22 19 16 13 10 7 4 32 26 20 14 8 2 1 5 11 17 23 29 3 6 9 12 15 18 21 24 27 30];
intanInds = intan(arrayfun(@(x) find(neuro==x), connectorInds));

% get desired site locations map
siteLocations = nan(32, 2);
siteLocations(4:3:31, :) = cat(2, ones(10,1)*-18, fliplr(12.5 : 25 : (12.5+25*9))'); % xy values for left column
siteLocations([1 2:3:32], :) = cat(2, zeros(12,1), fliplr(0 : 25 : (25*11))'); % xy values for right column
siteLocations(3:3:30, :) = cat(2, ones(10,1)*18, fliplr(12.5 : 25 : (12.5+25*9))'); % xy values for right column

% visualize final mapping results
% close all; figure;
% scatter(siteLocations(:,1), siteLocations(:,2))
% for i = 1:size(siteLocations,1)
%     text(siteLocations(i,1), siteLocations(i,2), num2str(i))
% end

% apply channel mapping to site locations
siteLocationsRemapped = nan(size(siteLocations));
for i = 1:32
    shankLocation = siteLocations(kiloSortInds(i),:);
    siteLocationsRemapped(intanInds(i),:) = shankLocation;
end



Nchannels = 32;
connected = true(Nchannels, 1);
chanMap   = 1:Nchannels;
chanMap0ind = chanMap - 1;
xcoords   = siteLocationsRemapped(:,1);
ycoords   = siteLocationsRemapped(:,2);
kcoords   = ones(Nchannels,1); % grouping of channels (i.e. tetrode groups)

fs = 30000; % sampling frequency
save('Z:\RAW\obstacleData\ephys\other\channelMaps\poly32intan.mat', ...
    'chanMap','connected', 'xcoords', 'ycoords', 'kcoords', 'chanMap0ind', 'fs')
