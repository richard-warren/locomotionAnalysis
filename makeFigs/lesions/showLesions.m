
% settings
micronToPixelRatio = 13.48;
sectionThickness = .1; % mm

folder = 'C:\Users\rick\Desktop\MTC3set';
files = dir(fullfile(folder, '*.roi'));
files = {files.name};

%%


aps = fliplr(0:sectionThickness:sectionThickness*length(files));
xys = nan(length(files)*2, 2);

ind = 1;
for i = 1:length(files)
    roi = ReadImageJROI(fullfile(folder, files{i}));
    xMin = min(roi.mnCoordinates(:,1));
    xMax = max(roi.mnCoordinates(:,1));
    
    xMin = xMin*micronToPixelRatio/1000;
    xMax = xMax*micronToPixelRatio/1000;
    
    xys(ind,:) = [xMin aps(i)];
    xys(length(files)*2-i+1,:) = [xMax aps(i)];
    ind = ind+1;
end

close all; figure; plot(xys(:,1), xys(:,2))
daspect([1 1 1])

%%

mat = roi.mnCoordinates;
mat = mat*micronToPixelRatio/1000;
close all; figure; plot(mat(:,1), -mat(:,2)); daspect([1 1 1])





