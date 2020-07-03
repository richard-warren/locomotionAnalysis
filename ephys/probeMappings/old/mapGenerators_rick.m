%% GENERATE MAPS FOR KILOSORT (I THINK THIS IS FOR VISUALIZATION PURPOSES ONLY)




%% NeuronexusA4x16-Poly2, intan

kiloSortInds = [6 8 10 12 14 16 4 2 1 3 15 13 11 9 7 5 22 24 26 28 30 32 20 18 17 19 31 29 27 25 23 21 38 40 42 44 46 48 36 34 33 35 47 45 43 41 39 37 54 56 58 60 62 64 52 50 49 51 63 61 59 57 55 53];
intanInds = [47 43 42 39 38 37 45 34 41 32 36 49 35 51 33 53 48 55 50 57 52 60 54 62 56 58 63 61 59 44 46 40 22 16 18 5 3 1 4 6 0 8 2 10 7 12 9 14 11 31 13 29 15 26 30 23 28 19 27 24 25 20 21 17] + 1;
channelNum_OpenEphys = [42 35 33 46 54 48 34 44 52 43 36 40 50 39 37 38 57 63 59 55 41 49 47 56 45 51 60 58 62 53 64 61 1 7 9 5 15 23 10 17 13 19 8 6 11 4 3 2 29 24 20 31 18 12 22 32 21 14 26 30 25 16 28 27];

% get desired site locations map
shankSeparation = 1000;
siteLocations = nan(16, 2); % start with one shank, then replicate after
siteLocations(2:2:16,:) = cat(2, zeros(8,1), (0:46:46*7)'+23); % first column of leftmost shank
siteLocations(1:2:15,:) = cat(2, ones(8,1)*30, (0:46:46*7)'); % second column of leftmost shank

siteLocations = repmat(siteLocations,4,1);
siteLocations(:,1) = siteLocations(:,1) + repelem(0:3,16)'.*shankSeparation;


siteLocations(17:32, 2) = siteLocations(17:32, 2) - 370; 
siteLocations(33:48, 2) = siteLocations(33:48, 2) - 370*2;
siteLocations(49:64, 2) = siteLocations(49:64, 2) - 370*3;






% apply channel mapping to site locations
siteLocationsRemapped = nan(size(siteLocations));
for i = 1:64
    shankLocation = siteLocations(kiloSortInds(i),:);
    siteLocationsRemapped(intanInds(i),:) = shankLocation;
end


Nchannels = 64;
connected = true(Nchannels, 1);
chanMap   = 1:Nchannels;
chanMap0ind = chanMap - 1;
xcoords   = siteLocationsRemapped(:,1);
ycoords   = siteLocationsRemapped(:,2);
kcoords_1   = repmat(1, 16, 1);
kcoords_2   = repmat(2, 16, 1);
kcoords_3   = repmat(3, 16, 1);
kcoords_4   = repmat(4, 16, 1);
kcoords = [kcoords_1; kcoords_2; kcoords_3; kcoords_4]; % grouping of channels (i.e. tetrode groups)
% grouping of channels (i.e. tetrode groups)
fs = 30000; % sampling frequency


% probe BDFD // original mapping
save('Y:\obstacleData\ephys\channelMaps\kilosort\BDFD.mat', ...
    'chanMap','connected', 'xcoords', 'ycoords', 'kcoords', 'chanMap0ind', 'fs')
% Z:\obstacleData\ephys\channelMaps\kilosort\forQualityMetricsPlot

% % probe BDFD2 // after right most shank broken
% connected([27 16 30 14 32 12 31 24 28 25 26 21 22 18 20 29]) = false;
% save('Z:\obstacleData\ephys\channelMaps\kilosort\BDFD2.mat', ...
%     'chanMap','connected', 'xcoords', 'ycoords', 'kcoords', 'chanMap0ind', 'channelNum_OpenEphys', 'fs')


% visualize final mapping results

scale = 1;

close all; figure('Color', 'white', 'Position', get(0,'ScreenSize'));
scatter(siteLocations(:,1)*scale, siteLocations(:,2)*scale)
for i = 1:size(siteLocations,1)
    text(siteLocations(i,1)*scale, siteLocations(i,2)*scale, [num2str(i) ' (' num2str(intanInds(kiloSortInds==i)-1) ')'])
end
daspect([1 1 1])



% probe BDFD // after right most shank broken
connected([27 16 30 14 32 12 31 24 28 25 26 21 22 18 20 29]) = false;
save('Y:\obstacleData\ephys\channelMaps\kilosort\BDFD2.mat', ...
    'chanMap','connected', 'xcoords', 'ycoords', 'kcoords', 'chanMap0ind', 'fs')


%% NeuronexusA1x32Poly3, intan

channelNum_OpenEphys = [1 17 8 24 16 2 29 32 7 26 3 15 21 19 11 23 14 12 28 30 6 18 9 13 22 25 5 27 10 4 31 20];


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

% % probe C6CE
% save('Z:\obstacleData\ephys\channelMaps\kilosort\C6CE.mat', ...
%     'chanMap', 'channelNum_OpenEphys', 'connected', 'xcoords', 'ycoords', 'kcoords', 'chanMap0ind', 'fs')

% probe D55F // nine sites defective?
connected([18 28 14 9 10 2 7 15 11]) = false;
save('Z:\obstacleData\ephys\channelMaps\kilosort\D55F.mat', ...
    'chanMap', 'channelNum_OpenEphys', 'connected', 'xcoords', 'ycoords', 'kcoords', 'chanMap0ind', 'fs')

%% NeuroCambridge ASSY-H2 (32 channel per shank*2)

% channelNum_OpenEphys is the order of openephys output .continuous files reflecting the
% real physical location of the sites on the probe. First 32 is shankA and
% 33-64 is shankB (refer to the cambridge neurotech map instruction)..
channelNum_OpenEphys = [28 26 24 19 21 15 32 30 18 29 31 22 20 23 25 27 1 3 6 8 10 12 14 16 17 13 11 9 7 5 4 2 63 61 60 58 56 54 52 50 47 51 53 55 57 59 62 64 38 40 42 45 43 49 34 36 48 35 33 44 46 41 39 37];
                    
channelNum_OpenEphys_shank1 = [28 26 24 19 21 15 32 30 18 29 31 22 20 23 25 27 1 3 6 8 10 12 14 16 17 13 11 9 7 5 4 2 ];
channelNum_OpenEphys_shank2 = [63 61 60 58 56 54 52 50 47 51 53 55 57 59 62 64 38 40 42 45 43 49 34 36 48 35 33 44 46 41 39 37];

y = [775:-25:0]'; % y coordinates for both shanks.
x_shank1 = repmat(1, 32, 1); % x coordinates for shankA.
x_shank2 = repmat(100, 32, 1); % x coordinates for shankB.
x = [x_shank1 x_shank2]; % x coordinates for both shanks.
sitaLocation_shank1 = [x_shank1 y];
siteLocation_shank2 = [x_shank2 y];
siteLocation = [sitaLocation_shank1; siteLocation_shank2]; % xy coordinates for both shanks.

% plot the probe channel map
figure;
scatter(siteLocation(:,1), siteLocation(:,2))
for i = 1:size(siteLocation,1)
    text(siteLocation(i,1), siteLocation(i,2), [num2str(i) ' (' num2str(channelNum_OpenEphys(i)-1) ')']) 
end

% save the channel map for kilosort
kccords_shank1 = ones(32, 1); 
kccords_shank2 = repmat(2, 32, 1); 

Nchannels = 64; 
connected = true(Nchannels, 1); 
chanMap = 1:Nchannels; 
chanMap0ind = chanMap - 1; 

for i = 1:64
    temp = find(channelNum_OpenEphys == i);
    xcoords(i) = siteLocation(temp, 1);
    ycoords(i) = siteLocation(temp, 2);
end

xcoords = xcoords';
ycoords = ycoords';

kcoords   = [kccords_shank1; kccords_shank2]; % grouping of channels (i.e. tetrode groups)
fs = 30000; % sampling frequency

% probe ASSY77
save('Z:\obstacleData\ephys\channelMaps\kilosort\forQualityMetricsPlot\ASSY77.mat', ...
    'chanMap','channelNum_OpenEphys','connected', 'xcoords', 'ycoords', 'kcoords', 'chanMap0ind', 'fs')





%% NeuroNexus A1x32-Ploy2 C585

channelNum_OpenEphys = [8, 24, 2, 29, 7, 26, 15, 21, 11, 23, 12, 28, 6, 18, 13, 22, 5, 27, 4, 31, 10, 20, 9, 25, 14, 30, 3, 19, 16, 32, 1, 17];

y1 = [850:-50:100]';
y2 = [825:-50:75]';
x2 = repmat(0, 16, 1);
x1 = repmat(43.3, 16, 1);

x = [];
y = [];
for i = 1:16
   x = [x; x1(i); x2(i)];
   y = [y; y1(i); y2(i)];
    
end
siteLocation = [x, y];
clear x1 x2 y1 y2


Nchannels = 32; 
connected = true(Nchannels, 1); 
chanMap = 1:Nchannels; 
chanMap0ind = chanMap - 1; 

for i = 1:32
    temp = find(channelNum_OpenEphys == i);
    xcoords(i) = siteLocation(temp, 1);
    ycoords(i) = siteLocation(temp, 2);
end

xcoords = xcoords';
ycoords = ycoords';
kcoords = ones(Nchannels,1);

fs = 30000; % sampling frequency

% plot the probe channel map
figure;
scatter(siteLocation(:,1), siteLocation(:,2))
for i = 1:size(siteLocation,1)
    text(siteLocation(i,1), siteLocation(i,2), [num2str(i) ' (' num2str(channelNum_OpenEphys(i)-1) ')']) 
end


% probe C858
save('Z:\obstacleData\ephys\channelMaps\kilosort\C858.mat', ...
    'chanMap','channelNum_OpenEphys','connected', 'xcoords', 'ycoords', 'kcoords', 'chanMap0ind', 'fs')

