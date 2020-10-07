ccf = loadCCF();
vol = ccf.coarseLabels==2;
vol = vol(1:10:end, 1:10:end, 1:10:end);  % downsample

%%
[x,y,z] = ind2sub(size(vol), find(vol));
xyz = boundary(x,y,z);
[X,Y] = meshgrid(xyz(:,1), xyz(:,2));
close all; figure;
% contour3(X,Y,Z)

%%

[Xs,Ys,Zs] = sphere(5);
close all; figure;
contour3(Xs, Ys, Zs)

%%
close all
lines = 25;
close all; figure('color', 'white', 'position', [680.00 141.00 1165.00 837.00]);
vol = double(ccf.coarseLabels(1:1:end, 1:1:end, 1:1:end));
vol = vol>0;

% tic
% for i = 1:size(vol,1)
% %     vol(i,:,:) = imfill(boundarymask(squeeze(vol(i,:,:))), 'holes');
%     [x,y] = ind2sub(size(vol,[2 3]), find(squeeze(vol(i,:,:))));
%     k = boundary(x,y);
%     vol(i,:,:) = poly2mask(x(k), y(k), size(vol,2), size(vol,3));
% end
% toc

vol = imclose(vol, strel('sphere', 10));
vol = imfill(vol, 'holes');
% vol = imfill(bwperim(vol), 'holes');
vol = smooth3(vol, 'box', 7);
vol = vol>.5;

mlslice = linspace(1,size(vol,3),lines);
apslice = linspace(1,size(vol,1),lines);
dvslice = []; %linspace(1,size(vol,2),lines);
thang = contourslice(vol, dvslice, apslice, mlslice, [.5 .5], 'cubic');
colormap gray
set(gca, 'zdir', 'reverse')
axis equal;
view(-45, 30)
figure; montage(vol)


%%

vol = double(ccf.coarseLabels(1:2:end, 1:2:end, 1:2:end));
vol = vol>0;

boundaries = boundary(vol);
mask = poly2mark();
% poly









