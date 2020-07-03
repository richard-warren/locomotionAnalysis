function results = ifInDCN(units, DCN, opts)


s.buffer = 150;
s.fileName = fullfile(getenv('OBSDATADIR'), 'ephys', 'ephysHistoData', 'ephysHistoData.mat');

% reassign settings contained in opts
if exist('opts', 'var'); for i = 1:2:length(opts); s.(opts{i}) = opts{i+1}; end; end


DCN_APCoords = unique(DCN(:, 2));
DCN_MLCoords = unique(DCN(:, 1));
DCN_DVCoords = unique(DCN(:, 3));

results = nan(size(units, 1), 1);
for i = 1:size(units, 1)
    
    % find the adjacent brain section
    if units(i, 2) > (DCN_APCoords(end) + s.buffer) || units(i, 2) < (DCN_APCoords(1) - s.buffer)
        results(i, 1) = false;
        continue
    else
        adjacentAPCoord = DCN_APCoords(knnsearch(DCN_APCoords, units(i, 2)), 1);
        
        %determine the DV and ML range
        section = DCN(find(DCN(:, 2) == adjacentAPCoord), :);
        section_MLCoords = unique(section(:, 1));
        section_DVCoords = unique(section(:, 3));
        x = section_MLCoords(knnsearch(section_MLCoords, units(i, 1)), 1);
        z = section_DVCoords(knnsearch(section_DVCoords, units(i, 3)), 1);
        
        
        inds = find(section(:, 3) == z);
        temp = sort(section(inds, 1));
        xmin = temp(1);
        xmax = temp(end);
        
        inds = find(section(:, 1) == x);
        temp = sort(section(inds, 3));
        zmin = temp(1);
        zmax = temp(end);
        
        % determine if unit is within the range 
        if units(i, 1) <= (xmax+s.buffer) && units(i, 1) >= (xmin-s.buffer) && units(i, 3) <= (zmax+s.buffer) && units(i, 3) >= (zmin-s.buffer)
            results(i, 1) = true;
        else
            results(i, 1) = false;
        end
    end    
    
    
end








 





end