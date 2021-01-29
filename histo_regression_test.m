%% fit linear model relating manipulator and histological coordinate systems

%% load data

load('C:\Users\richa\Desktop\manipuPoints.mat')
load('C:\Users\richa\Desktop\histoPoints.mat')

%% make matrices of of manipulator and histological coordinates

dii_bins  = ~contains(manipuPoints.ManipuProbeNames, 'noDi');
[~, inds] = ismember(manipuPoints{dii_bins,'ManipuProbeNames'}, histoPoints.HistoProbeNames);  % indices for histoPoints that match the dii traces in manipPoints

manip       = manipuPoints.ManipuCoords(dii_bins, :);   % manipulator coordinates for dii traces (nx2)
histo       = histoPoints.HistoCoords(inds, :);         % histo coordinates for dii traces
manip_nodii = manipuPoints.ManipuCoords(~dii_bins, :);  % manipulator coordinates for no dii traces

n_dii   = size(manip, 1);        % number of dii tracks
n_nodii = size(manip_nodii, 1);  % number of no dii traces

%% we want to find a matrix that maps manipulator -> histo coordinates!
%  first lets pack everything into nice matrices / vectors

M       = [manip, ones(n_dii,1)];          % pad the manipulator coordinates with ones - this will allow us to fit a bias term
M_nodii = [manip_nodii, ones(n_nodii,1)];  % do the same for the no dii traces
h = histo;

%% now we want to solve M X = h, where M is manipulator coordinates, h is histo coordinates, 
%  and X is a matrix that defines the linear transformation

X = M \ h;  % fancy one line matlab syntax! this solves the linear equation MX = h for the matrix X

%% now we can generate predictions for the histo coordinates, both for the 
%  'training' data (dii traces) and the 'test' data (no dii traces)
histo_predicted = M * X;
nodii_predicted = M_nodii * X;

%% plot!

figure('color', 'white', 'position', [512.00 815.00 678.00 250.00])
colors = lines(n_dii);

subplot(1,2,1); hold on
title('maniupulator');
scatter(manip(:,1), manip(:,2), 100, colors)
scatter(manip_nodii(:,1), manip_nodii(:,2), 100, 'black', 'filled')
legend('dii', 'no dii', 'Location', 'best')

subplot(1,2,2); hold on
title('histology');
scatter(histo(:,1), histo(:,2), 100, colors)
scatter(histo_predicted(:,1), histo_predicted(:,2), 100, colors, '+')
scatter(nodii_predicted(:,1), nodii_predicted(:,2), 100, 'black', 'filled')
legend('dii (actual)', 'dii (predicted)', 'no dii')








