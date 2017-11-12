function filled = fillShortMissing(data, maxGapToFill, method)

    % like fillmissing (matlab function), but only fills short gaps in data, with max gap width of maxGapTofill
    
    % orient data horizontally
    if size(data,1)>size(data,2)
        flipped = true;
        data = data';
    else
        flipped = false;
    end
    
    % fill all values
    filled = fillmissing(data, method, 'EndValues', 'none');

    % find indices where nan epohcs begin and end
    nanStartInds = find(diff(isnan(data))==1)+1;
    nanStopInds = find(diff(isnan(data))==-1);

    % if data start or end with nan, add this to start or end index lists
    if isnan(data(1)); nanStartInds = [1 nanStartInds]; end
    if isnan(data(end)); nanStopInds = [nanStopInds length(data)]; end
    
    % find indices of nan periods that are greater than maxGapToFill
    nanLengths = nanStopInds - nanStartInds;
    longNanInds = find(nanLengths > maxGapToFill);

    % restore NaN values to all periods greater than maxGapToFill
    for j = longNanInds
        filled(nanStartInds(j):nanStopInds(j)) = nan;
    end
    
    % flip data back to original orientation
    if flipped; filled = filled'; end

end