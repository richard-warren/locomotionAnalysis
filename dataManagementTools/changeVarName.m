function changeVarName(originalName, newName)

% change name of variable saved in run.mat from originalName to newName

dataDir = 'C:\Users\Rick\Google Drive\columbia\obstacleData\sessions\';
dataFolders = dir(dataDir);
dataFolders = dataFolders(3:end); % remove current and parent directory entries
dataFolders = dataFolders([dataFolders.isdir]); % keep only folders

for i = 1:length(dataFolders)
    
    % load data
    file = [dataDir dataFolders(i).name '\run.mat'];
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