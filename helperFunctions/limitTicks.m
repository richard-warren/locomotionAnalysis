function limitticks()
% removes all but two tick labels for y axis of current figure
% todo: add x or y limits only option

pause(.001)

xticks = limitticks_(get(gca, 'XTickLabel'));
set(gca, 'XTickLabel', xticks)

yticks = limitticks_(get(gca, 'YTickLabel'));
set(gca, 'YTickLabel', yticks)


function ticks = limitticks_(ticks)
    labelBins = ~cellfun(@isempty, ticks);  % bins where there are tick labels
    yticksSub = ticks(labelBins);

    if sum(labelBins)>2

        yticksNum = cellfun(@str2num, yticksSub);  % convert to nums

        % if zero is present, use 0 and whichever label is furthest from zero
        % axis is symmetric around zero, keep biggest label on both sides of zero
        if any(yticksNum==0)
            absMax = max(abs(yticksNum));

            % set all but zero and max ind to empty arrays
            inds = 1:length(yticksSub);
            inds(ismember(abs(yticksNum), [0 absMax])) = [];
            yticksSub(inds) = {''};

        % otherwise just use first and last tick
        else
            yticksSub(2:end-1) = {''};
        end
    end

    ticks(labelBins) = yticksSub; 
end

end