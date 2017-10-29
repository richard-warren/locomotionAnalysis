function pimpFig
    
    % make window big
    gcf
    set(gcf, 'units', 'normalized', 'position', [0 0 1 1], 'color', [1 1 1]) 
    
    % adjust axis settings
    axes = findobj(gcf,'type','axes');
    for a = axes
        set(a, 'box', 'off');
    end
end