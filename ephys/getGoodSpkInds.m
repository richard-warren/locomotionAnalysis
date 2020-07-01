function [spkInds, unit_ids] = getGoodSpkInds(session)

% given session name, returns cell array where each entry is vector of
% spike times for units marked as 'good' in kilosort


% initializations
ephysFolder = dir(fullfile(getenv('OBSDATADIR'), 'sessions', session, 'ephys_*'));
ephysFolder = fullfile(ephysFolder.folder, ephysFolder.name);
addpath(fullfile(getenv('GITDIR'), 'npy-matlab'))
addpath(fullfile(getenv('GITDIR'), 'analysis-tools'))
allSpkInds = readNPY(fullfile(ephysFolder, 'spike_times.npy'));
clusters = readNPY(fullfile(ephysFolder, 'spike_clusters.npy'));


% find good units
if exist(fullfile(ephysFolder, 'cluster_groups.csv'), 'file')  % old kilosort1 format
    clusterGroups = readtable(fullfile(ephysFolder, 'cluster_groups.csv'));
    unit_ids = clusterGroups.cluster_id(strcmp(clusterGroups.group, 'good'));
    
elseif exist(fullfile(ephysFolder, 'cluster_group.tsv'), 'file')  % new kilosort1 format
    clusterGroups = tdfread(fullfile(ephysFolder, 'cluster_group.tsv'));
    unit_ids = clusterGroups.cluster_id(all(clusterGroups.group=='good',2));
end


% get spike times for individual units
spkInds = cell(1,length(unit_ids));
for i = 1:length(unit_ids)
    spkInds{i} = allSpkInds(clusters==unit_ids(i));
end








