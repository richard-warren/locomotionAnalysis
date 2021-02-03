function limitTicks()
% removes all but two tick labels for y axis of current figure

pause(.001)
yticks = get(gca, 'YTickLabel');
labelBins = ~cellfun(@isempty, yticks);  % bins where there are tick labels
yticksSub = yticks(labelBins);

if sum(labelBins)>2
    
    yticksNum = cellfun(@str2num, yticksSub);  % convert to nums
    
    % if zero is present, use 0 and whichever label is furthest from zero
    zeroInd = find(yticksNum==0);
    
    if ~isempty(zeroInd)
        [~, maxInd] = max(abs(yticksNum));
        
        % set all but zero and max ind to empty arrays
        inds = 1:length(yticksSub);
        inds([zeroInd maxInd]) = [];
        yticksSub(inds) = {''};
    
    % otherwise just use first and last tick
    else
        yticksSub(2:end-1) = {''};
    end
end

yticks(labelBins) = yticksSub;
set(gca, 'YTickLabel', yticks)