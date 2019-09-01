% for vgat opto experiments, i shine a light on the skull and activate vgat
% neurons transcranially // the diameter of the light spot on the skull is
% related to laser path cable diam, NA, and distance to skull, among
% other potential factors // svaboda paper reports 400um diameter, when
% diameter is measured as 4 * (standard deviation) // this code measures
% light spot diam similarly so i can compare to guo paper...

% settings
% clear all; close all
imgFile = 'Y:\obstacleData\other\lightDiamPics\2mm tests\attempt2\4.bmp';
pixPerMm = 318;

% read image
img = imread(imgFile); img = img(:,:,3);  % blue light, so take only blue channel

% get light intensity in brightest row
[~, maxRow] = max(mean(img,2));  % find brightest row, on which subsequent analyses will be performed
intensity = double(img(maxRow,:)); % intensity for that row
intensity = intensity / sum(intensity);

% find mean (first moment) and std (sqrt of second moment)
firstMoment = sum((1:length(intensity)).*intensity);
x = ((1:size(img,2)) - firstMoment) / pixPerMm;

secondMoment = sum((1:length(intensity)).^2.*intensity);
stDev = sqrt(secondMoment - firstMoment^2);




% plot image of light spot
figure('position', [1962 68 500 900]);
subplot(2,1,1);
imagesc(img); hold on
scatter(firstMoment, maxRow, 'red', '+')

% plt row intensity
subplot(2,1,2);
plot(x, intensity); hold on
line([0 0], get(gca, 'YLim'), 'color', 'red')
line(([-stDev stDev])/pixPerMm, [0 0], 'color', 'red', 'linewidth', 3)
set(gca, 'xlim', [x(1) x(end)])


lightDiam = stDev*4/pixPerMm;  % definition of light diam from guo, svaboda...
fprintf('estimated light diameter: %.2f\n', lightDiam)
