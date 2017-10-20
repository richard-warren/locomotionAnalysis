function changeVarName(originalName, newName, matFile)

% goes through all folders in dataDir, and in matFile changes the name of variable originalName to newName
%
% input       originalName:   original name of variable contained within matFile
%             newName:        name to which original name will be changed
%             matFile:        name of .mat file containing the variable to be changed (e.g. 'run.mat')


dataDir = 'C:\Users\Rick\Google Drive\columbia\obstacleData\sessions\';
dataFolders = dir(dataDir);
dataFolders = dataFolders(3:end); % remove current and parent directory entries
dataFolders = dataFolders([dataFolders.isdir]); % keep only folders

for i = 1:length(dataFolders)
    
    % load data
    file = [dataDir dataFolders(i).name '\' matFile];
    data = load(file);
    
    if any(strcmp(fieldnames(data), originalName))
        
        fprintf(['changing ''' originalName ''' to ''' newName ''' in %s\n'], dataFolders(i).name)
        
        % create copy of original field
        eval(['data.' newName ' = data.' originalName ';'])

        % delete original field
        data = eval(['rmfield(data, ''' originalName ''')']);

        % save file
        save(file, '-struct', 'data');
    else
        fprintf(['''' originalName ''' not found in %s\n'], dataFolders(i).name)
    end
    
end