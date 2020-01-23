% PSTHs center around swing start and stance start timepoints, 
    % for each paw
    for paw = 1:4
        % plot PSTH center around swing start point
        plotInd = plotInd + 1; subplot(s.rows, s.cols, plotInd);
        
        swingStartInds = find(diff(~stanceBins(:,paw))==1);
        swingStartTimes = frameTimeStamps(swingStartInds);        
        validBins = find(swingStartTimes > minTime & swingStartTimes < maxTime);
        validSwingStartTimes = swingStartTimes(validBins);
        
        swingEndInds = find(diff(stanceBins(:,paw))== 1);
        swingEndTimes = frameTimeStamps(swingEndInds);              
        validBins = find(swingEndTimes > minTime & swingEndTimes < maxTime);
        validSwingEndTimes = swingEndTimes(validBins);
        
        % quality check
        if length(validSwingStartTimes) > length(validSwingEndTimes)
            if validSwingStartTimes(end) > validSwingEndTimes(end)
                validSwingStartTimes = validSwingStartTimes(1:end-1);
            else
                validSwingStartTimes = validSwingStartTimes(2:end);
            end
            
        elseif length(validSwingStartTimes) < length(validSwingEndTimes)
            if validSwingEndTimes(1) < validSwingStartTimes(1)
                validSwingEndTimes = validSwingEndTimes(2:end);
            else
                validSwingEndTimes = validSwingEndTimes(1:end-1);
            end
        end
        
        substracted = validSwingEndTimes - validSwingStartTimes;
        if any(substracted < 0)
            fprintf('WARNING: swing start and end timepoints mismatch!');
        end
        
        % discard steps which swing phases are unreasonably long
        if any(substracted > 0.5)
            validSwingEndTimes(find(substracted > 0.5)) = [];
            validSwingStartTimes(find(substracted > 0.5)) = [];
        end
        
        % plot PSTH center around swing start point                    
        plotPSTH2(session, cellNum, validSwingStartTimes, {'xLims', [-0.05 0.05], 'colors', s.pawColors(paw,:)});
        xlabel(sprintf('%s: swing start', s.pawNames{paw}))
        
        % plot PSTH center around stance start point
        plotInd = plotInd + 1; subplot(s.rows, s.cols, plotInd);               
        plotPSTH2(session, cellNum, validSwingEndTimes, {'xLims', [-0.05 0.05], 'colors', s.pawColors(paw,:)});
        xlabel(sprintf('%s: stance start', s.pawNames{paw}))
    end