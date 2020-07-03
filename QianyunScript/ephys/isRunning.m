function [isRunning, selectedStartInds, selectedEndInds] = isRunning(frameVel, opts)

% settings
s.velThresh = 0.25;
s.frameFreq = 250; % Hz
s.durationThresh = 1.5; % second

% reassign settings contained in opts
if exist('opts', 'var'); for i = 1:2:length(opts); s.(opts{i}) = opts{i+1}; end; end

inds = frameVel > s.velThresh;
startInds = find(diff(inds) == 1) + 1;
endInds = find(diff(inds) == -1);

while(endInds(1) < startInds(1))
    endInds(1) = [];
end

while(startInds(end) > endInds(end))
    startInds(end) = [];
end

if length(startInds) ~= length(endInds)
    warning('UNEQUAL: length of startInds and endInds');
else
    durations = endInds - startInds;
    selected = find(durations > s.frameFreq*s.durationThresh);
    selectedStartInds = startInds(selected);
    selectedEndInds = endInds(selected);
    
    isRunning = nan(1, length(frameVel));
    isRunning(1, 1:selectedStartInds(1)) = 0;
    for i = 1:length(selectedStartInds)-1
        isRunning(1, selectedStartInds(i):selectedEndInds(i)) = 1;
        isRunning(1, selectedEndInds(i)+1:selectedStartInds(i+1)) = 0;
    end

    isRunning(1, selectedStartInds(end):selectedEndInds(end)) = 1;
    isRunning(1, selectedEndInds(end)+1:length(frameVel)) = 0;
  
end
        
        


end