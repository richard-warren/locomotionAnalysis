function [ids, names] = getAllDescenants(tree_json, parent)
% given the json file for the tree structure of the allen brain atlas,
% finds the names and ids of all children of 'parent' // use this to find
% all labels belonging to the cerebellum, for example

% 'tree_json' is downloaded from: https://community.brain-map.org/t/allen-mouse-ccf-accessing-and-using-related-data-and-tools/359

tree = jsondecode(fileread(tree_json));
tree = tree.msg;
parentStruct = getParentStruct(tree, parent);
[names, ids] = findDescendants(parentStruct);


function parentStruct = getParentStruct(struct, parent)
    % recursively find sub-struct within struct for desired parent struct
    
    for i = 1:length(struct)
        if strcmp(struct(i).name, parent)
            parentStruct = struct(i);
            break
        else
            % search down tree
            if isstruct(struct(i).children)
                parentStruct = getParentStruct(struct(i).children, parent);
                if isstruct(parentStruct); break; end
            else
                parentStruct = false;
            end
        end
    end
end

function [names, ids] = findDescendants(struct)
    % recursively find all descendenants of struct
    
    names = {struct.name};
    ids = [struct.id];
    
    for i = 1:length(struct)
        if isstruct(struct(i).children)
            [desNames, desIds] = findDescendants(struct(i).children);
            names = [names, desNames];
            ids = [ids, desIds];
        end
    end
end

end