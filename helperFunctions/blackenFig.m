function blackenFig


% settings
black = repmat(0, 1,3);
white = repmat(.85, 1,3);

colorSettings.figure = {{'Color', black}};
colorSettings.axes = {{'Color', black}, {'XColor', 'inverted'}, {'YColor', 'inverted'}};
colorSettings.line = {{'Color', 'inverted'}};
% colorSettings.line = {{'Color', white}};
colorSettings.scatter = {{'MarkerFaceColor', 'inverted'}};
colorSettings.text = {{'Color', 'inverted'}};
colorSettings.patch = {{'FaceColor', 'inverted'}};
colorSettings.rectangle = {{'FaceColor', 'inverted'}};
colorSettings.legend = {{'color', 'inverted'}, {'textcolor', 'inverted'}};


% initializations
obs = findall(gcf); % get all graphics objects
objectsWithSettings = fieldnames(colorSettings);



for i = 1:length(obs)
    
    if any(strcmp(obs(i).Type, objectsWithSettings)) % see if this object has settings listed above
        
        for j = 1:length(colorSettings.(obs(i).Type)) % iterate through settings for this object
            
            settingName = colorSettings.(obs(i).Type){j}{1};
            targetSetting = colorSettings.(obs(i).Type){j}{2};
            
                
            if strcmp(targetSetting, 'inverted')

                settingValue = get(obs(i), settingName);

                % this is the ob color for scatter data
                if strcmp(settingValue, 'flat') || strcmp(settingValue, 'none')
                    for k = 1:size(obs(i).CData,1)
                        settingValue = obs(i).CData(k,:);
                        if all(diff(settingValue)==0); obs(i).CData(k,:) = 1-settingValue; end % only invert if original color was gray
                    end
                % otherwise invert the color normally
                else
                    if all(diff(settingValue)==0); set(obs(i), settingName, 1-settingValue); end % only invert if original color was gray
                end

            % set to target color    
            else
                set(obs(i), settingName, targetSetting);
            end
            

        end
    end
end