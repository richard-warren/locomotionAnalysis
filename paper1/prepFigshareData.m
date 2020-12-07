%% prepare datasets for upload to open-source figshare data repository


%% load experiment data

% settings
filename = 'E:\google_drive\columbia\papers\hurdles_paper1\elife_submission\revision\obstacle_data.mat';

fprintf('loading... '); load(fullfile(getenv('SSD'), 'paper1', 'baseline_data.mat'), 'data'); disp('baseline data loaded!')
data = data.data;  % remove top level of hierarchy
save(filename, 'data')