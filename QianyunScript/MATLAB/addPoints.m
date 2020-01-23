function addPoints(model, d_PC, d_channels, d_start, d_end)



avg = model.avg;
dirVect = model.dirVect;

if exist('crossPoint_PC', 'var')
    crossPoint = (d_PC)*dirVect + avg;
    hold on
    plot3(crossPoint(:, 1), crossPoint(:, 3), crossPoint(:, 2), '.r', 'MarkerSize', 30)
end

if exist('d_start', 'var')
    startPoint = (d_start)*dirVect + avg;
    hold on
    plot3(startPoint(:, 1), startPoint(:, 3), startPoint(:, 2), '.m', 'MarkerSize', 20)
    
end


if exist('d_end', 'var')
    endPoint = (d_end)*dirVect + avg;
    hold on
    plot3(endPoint(:, 1), endPoint(:, 3), endPoint(:, 2), '.m', 'MarkerSize', 20)
    
end

if exist('d_channels', 'var')
    for i = 1:length(d_channels)
         channelPoint = d_channels(i)*dirVect + avg;
         hold on
         plot3(channelPoint(:, 1), channelPoint(:, 3), channelPoint(:, 2), '.c', 'MarkerSize', 10)        
    end
end




end
