% perform analyses that don't require access to engram

sessions = getEphysSessions();
% sessions = sessions(1:33);  % temp
% sessions = {'181020_001'};

overwrite = true;


tic
parfor i = 1:length(sessions)
    folder = fullfile(getenv('SSD'), 'paper2', 'modelling');
    try
        % neural responses
        filename = fullfile(folder, 'responses', [sessions{i} '_responses.mat']);
        if overwrite || ~exist(filename, 'file')
            getNeuralResponses(sessions{i})
        end
        
        % plot neural responses
        plotNeuralResponses(sessions{i}, 'visible', false, 'showImportance', false)

        % design matrices
%         filename = fullfile(folder, 'designMatrices', [sessions{i} '_designMatrix.mat']);
%         if overwrite || ~exist(filename, 'file')
%             makeDesignMatrix(sessions{i}, 'timeDegrees', 3, 'outputFileName', filename);
%         end
        
    catch exception
        fprintf('%s: PROBLEM! -> %s\n', sessions{i}, exception.identifier)
    end
end
fprintf('finished in %.1f minutes\n', toc/60)


%% histo

load(fullfile(getenv('OBSDATADIR'), 'histology', '0_ephysHistoData', 'ephysHistoTable.mat'), 'ephysHistoTable')
mice = unique(ephysHistoTable.mouseID);

overwrite = false;

for i = 1:length(mice)
    filename = ['E:\lab_files\paper2\histo\registration\' mice{i} '_registration.mat'];
    if ~exist(filename, 'file') || overwrite
        registerBrain(mice{i});
    end
end

%% count total units in each nucleus

nucleus = cell(1, length(mice));
for i = 1:length(mice)
    try
        load(['E:\lab_files\paper2\histo\registration\' mice{i} '_registration.mat'], 'registration')
        nucleus{i} = registration.nucleus;
    catch expection
        fprintf('%s: problem -> %s\n', mice{i}, expection.identifier)
    end
end
nucleus = cat(1, nucleus{:});
n = length(nucleus);

fprintf('\n\nUNIT COUNTS\n')
fprintf('-----------\n')
fprintf('fastigial:    %3i/%i\n', sum(strcmp(nucleus, 'fastigial')), n)
fprintf('interpositus: %3i/%i\n', sum(strcmp(nucleus, 'interpositus')), n)
fprintf('dentate:      %3i/%i\n', sum(strcmp(nucleus, 'dentate')), n)
fprintf('other:        %3i/%i\n', sum(strcmp(nucleus, 'none')), n)
nucleiTotal = n-sum(strcmp(nucleus, 'none'));
fprintf('nuclei:       %3i/%i (%.1f%%)\n', nucleiTotal, n, nucleiTotal*100/n)



%% todo (automatically copy files from engram to local, skipping files that are already there...)


