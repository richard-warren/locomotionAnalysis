function [spkInds, unit_ids, bestChannels] = getGoodSpkInds_old(session)

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
channel_map = readNPY(fullfile(getenv('OBSDATADIR'), 'sessions', session, ephysFolder, 'channel_map.npy')) + 1; % add 1 to avoid zero indexing
bestChannels = nan(1,length(unit_ids));

for i = 1:length(unit_ids)
    
    % find templates used in cluster
    mostCommonTemplateInd = mode(spike_templates(clusters==unit_ids(i))) + 1; % add 1 to avoid zero indexing
    cluster_template = squeeze(templates(mostCommonTemplateInd,:,:));
    
    % average templates used in cluster, and find channel with max peak to peak waveform
%     mean_template = squeeze(mean(cluster_templates,1));
%     mean_template = squeeze(mean(cluster_templates,1));
    [~, ind] = max(peak2peak(cluster_template,1));
    bestChannels(i) = channel_map(ind);
end



% get spike times for individual units
spkInds = cell(1,length(unit_ids));
for i = 1:length(unit_ids)
    spkInds{i} = allSpkInds(clusters==unit_ids(i));
end