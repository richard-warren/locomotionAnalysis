function updateTextAndLine(text, line, eventTimes, textString)
    eventInd = find(eventTimes>=timeEpochs(trials(i),1)-s.voltageWindow & ...
                    eventTimes<=timeEpochs(trials(i),2),1,'first');
    if ~isempty(eventInd)
        eventTime = eventTimes(eventInd);
        set(text, 'Position', [eventTime+s.voltageWindow*.01 s.yLims(2)])
        set(line, 'XData', [eventTime eventTime])
    end
    if exist('textString', 'var'); set(text, 'String', textString); end
end