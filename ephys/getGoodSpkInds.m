function [spkInds, unit_ids, bestChannels] = getGoodSpkInds(session)

% given session name, returns cell array where each entry is vector of
% spike times for units marked as 'good' in kilosort

% get name of ephys folder
files = dir(fullfile(getenv('OBSDATADIR'), 'sessions', session));
ephysFolder = files([files.isdir] & contains({files.name}, 'ephys_')).name;


addpath(fullfile(getenv('GITDIR'), 'npy-matlab'))
addpath(fullfile(getenv('GITDIR'), 'analysis-tools'))
allSpkInds = readNPY(fullfile(getenv('OBSDATADIR'), 'sessions', session, ephysFolder, 'spike_times.npy'));
clusters = readNPY(fullfile(getenv('OBSDATADIR'), 'sessions', session, ephysFolder, 'spike_clusters.npy'));
clusterGroups = readtable(fullfile(getenv('OBSDATADIR'), 'sessions', session, ephysFolder, 'cluster_groups.csv'));
unit_ids = clusterGroups.cluster_id(strcmp(clusterGroups.group, 'good'));

% get spike times for individual units
spkInds = cell(1,length(unit_ids));
for i = 1:length(unit_ids)
    spkInds{i} = allSpkInds(clusters==unit_ids(i));
end

% get best channels
spike_templates = readNPY(fullfile(getenv('OBSDATADIR'), 'sessions', session, ephysFolder, 'spike_templates.npy'));
templates = readNPY(fullfile(getenv('OBSDATADIR'), 'sessions', session, ephysFolder, 'templates.npy'));
bestChannels = nan(1,length(unit_ids));
for i = 1:length(unit_ids)
    
    % find templates used in cluster
    cluster_templates_inds = unique(spike_templates(clusters==unit_ids(i))) + 1; % add 1 to avoid zero indexing
    cluster_templates = templates(cluster_templates_inds,:,:);
    
    % average templates used in cluster, and find channel with max peak to peak waveform
    mean_template = squeeze(mean(cluster_templates,1));
    [~, bestChannels(i)] = max(peak2peak(mean_template,1));
end
