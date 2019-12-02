function flat = flattenData_nonrecursive(data, vars)

% this is an ugly, non-recursive version of flattenData // i'm using this
% to check that both algorithms return the same values // it is actually
% faster when length(vars)<~45, but it is also just a really ugly algorithm
% with a lot of assumptions built in...


flat = struct();
ind = 1;
[mouseVars, sessionVars, trialVars, pawVars] = deal({});

mouseVars = vars(ismember(vars, fieldnames(data.data)));
sessionVars = vars(ismember(vars, fieldnames(data.data(1).sessions)));
trialVars = vars(ismember(vars, fieldnames(data.data(1).sessions(1).trials)));
pawVars = vars(ismember(vars, fieldnames(data.data(1).sessions(1).trials(1).paws)));

levels = {'mouse', 'session', 'trial', 'paw'};
deepestLevel = levels{find(~cellfun(@isempty, {mouseVars, sessionVars, trialVars, pawVars}), 1, 'last')};


for i = 1:length(data.data)
    if strcmp(deepestLevel, 'mouse')
        addData();
    else
        for j = 1:length(data.data(i).sessions)
            if strcmp(deepestLevel, 'session')
                addData();
            else
                for k = 1:length(data.data(i).sessions(j).trials)
                    if strcmp(deepestLevel, 'trial')
                        addData();
                    else
                        for l = 1:length(data.data(i).sessions(j).trials(k).paws)
                            addData();
                        end
                    end
                end
            end
        end
    end
end




function addData()
    try
        % store mouse vars
        for m = 1:length(mouseVars)
            flat(ind).(mouseVars{m}) = data.data(i).(mouseVars{m});
        end

        % store session vars
        for m = 1:length(sessionVars)
            flat(ind).(sessionVars{m}) = data.data(i).sessions(j).(sessionVars{m});
        end

        % store trial vars
        for m = 1:length(trialVars)
            flat(ind).(trialVars{m}) = data.data(i).sessions(j).trials(k).(trialVars{m});
        end

        % store paw vars
        for m = 1:length(pawVars)
            flat(ind).(pawVars{m}) = data.data(i).sessions(j).trials(k).paws(l).(pawVars{m});
        end

        ind = ind + 1;
    catch; end
end

end