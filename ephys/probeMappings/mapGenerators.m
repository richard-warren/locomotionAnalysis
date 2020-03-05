%% GENERATE MAPS FOR KILOSORT (I THINK THIS IS FOR VISUALIZATION PURPOSES ONLY)




%% NeuronexusA4x16-Poly2, intan

kiloSortInds = [6 8 10 12 14 16 4 2 1 3 15 13 11 9 7 5 22 24 26 28 30 32 20 18 17 19 31 29 27 25 23 21 38 40 42 44 46 48 36 34 33 35 47 45 43 41 39 37 54 56 58 60 62 64 52 50 49 51 63 61 59 57 55 53];
intanInds = [47 43 42 39 38 37 45 34 41 32 36 49 35 51 33 53 48 55 50 57 52 60 54 62 56 58 63 61 59 44 46 40 22 16 18 5 3 1 4 6 0 8 2 10 7 12 9 14 11 31 13 29 15 26 30 23 28 19 27 24 25 20 21 17] + 1;

% get desired site locations map
shankSeparation = 200;
siteLocations = nan(16, 2); % start with one shank, then replicate after
siteLocations(2:2:16,:) = cat(2, zeros(8,1), (0:46:46*7)'+23); % first column of leftmost shank
siteLocations(1:2:15,:) = cat(2, ones(8,1)*30, (0:46:46*7)'); % second column of leftmost shank

siteLocations = repmat(siteLocations,4,1);
siteLocations(:,1) = siteLocations(:,1) + repelem(1:4,16)'.*shankSeparation;


% apply channel mapping to site locations
siteLocationsRemapped = nan(size(siteLocations));
for i = 1:64
    shankLocation = siteLocations(kiloSortInds(i),:);
    siteLocationsRemapped(intanInds(i),:) = shankLocation;
end

% visualize final mapping results



Nchannels = 64;
connected = true(Nchannels, 1);
chanMap   = 1:Nchannels;
chanMap0ind = chanMap - 1;
xcoords   = siteLocationsRemapped(:,1);
ycoords   = siteLocationsRemapped(:,2);
kcoords   = ones(Nchannels,1); % grouping of channels (i.e. tetrode groups)
fs = 30000; % sampling frequency


% % probe BDFD // original mapping
% save('Y:\obstacleData\ephys\channelMaps\kilosort\BDFD.mat', ...
%     'chanMap','connected', 'xcoords', 'ycoords', 'kcoords', 'chanMap0ind', 'fs')

% probe BDFD // after right most shank broken
connected([27 16 30 14 32 12 31 24 28 25 26 21 22 18 20 29]) = false;
save('Z:\obstacleData\ephys\channelMaps\kilosort\BDFD2.mat', ...
    'chanMap','connected', 'xcoords', 'ycoords', 'kcoords', 'chanMap0ind', 'fs')



close all; figure('Color', 'white', 'Position', [134 119 1595 807]);
scatter(siteLocations(:,1)*0.1, siteLocations(:,2)*0.1)
for i = 1:size(siteLocations,1)
    text(siteLocations(i,1)*0.1, siteLocations(i,2)*0.1, [num2str(i) ' (' num2str(intanInds(kiloSortInds==i)-1) ')'])
end
daspect([1 1 1])



%% NeuronexusA1x32Poly3, intan

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


% apply channel mapping to site locations
siteLocationsRemapped = nan(size(siteLocations));
for i = 1:32
    shankLocation = siteLocations(kiloSortInds(i),:);
    siteLocationsRemapped(intanInds(i),:) = shankLocation;
end

% visualize final mapping results
close all; figure('Position', [134 119 400 807]);
scatter(siteLocations(:,1), siteLocations(:,2))
for i = 1:size(siteLocations,1)
    text(siteLocations(i,1), siteLocations(i,2), [num2str(i) ' (' num2str(intanInds(kiloSortInds==i)) ')'])
end


Nchannels = 32;
connected = true(Nchannels, 1);
chanMap   = 1:Nchannels;
chanMap0ind = chanMap - 1;
xcoords   = siteLocationsRemapped(:,1);
ycoords   = siteLocationsRemapped(:,2);
kcoords   = ones(Nchannels,1); % grouping of channels (i.e. tetrode groups)
fs = 30000; % sampling frequency

% probe C6CE
save('Y:\obstacleData\ephys\channelMaps\kilosort\C6CE.mat', ...
    'chanMap','connected', 'xcoords', 'ycoords', 'kcoords', 'chanMap0ind', 'fs')

% probe D55F // nine sites defective?
connected([18 28 14 9 10 2 7 15 11]) = false;
save('Y:\obstacleData\ephys\channelMaps\kilosort\D55F.mat', ...
    'chanMap','connected', 'xcoords', 'ycoords', 'kcoords', 'chanMap0ind', 'fs')



