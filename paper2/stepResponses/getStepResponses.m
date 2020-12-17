function data = getStepResponses()
    % creates table containing neural responses to all steps for all units in
    % all sessions

    % settings
    nx = 50;  % number of points in the x axis


    % inits
    sessions = getEphysSessions();
    data = table(cell(length(sessions),4), cell(length(sessions),4), cell(length(sessions),4), cell(length(sessions),1), ...
        'RowNames', sessions, 'VariableNames', {'stepData', 'responses', 'predicted', 'unit_ids'});
    outs = cell(1, length(sessions));  % stored outputs of parallelized fcn
    
    % compute sessions in parallel
    parfor i = 1:length(sessions)
        fprintf('(%2i/%i) getting step responses for session: %s\n', ...
            i, length(sessions), sessions{i})
        try
            outs{i} = getSesStepResponses(sessions{i}, nx);
        catch
            fprintf('PROBLEM WITH SESSION %s (%i)\n', sessions{i}, i)
        end
    end

    % unpack function outputs and stick into table
    for i = 1:length(sessions)
        data{i, 'stepData'} = outs{i}{1};
        data{i, 'responses'} = outs{i}{2};
        data{i, 'predicted'} = outs{i}{3};
        data{i, 'unit_ids'} = outs{i}(4);
    end
    
    disp('all done!')
end







function outs = getSesStepResponses(session, nx)
    % get step responses for single session

    fname = fullfile(getenv('SSD'), 'paper2', 'modelling', 'stepData', [session '_stepData.mat']);
    responses = cell(1,4);  % one response matrix per paw // each matrix is (unit X step X 'time')
    predicted = cell(1,4);
    stepData = cell(1,4);
    unit_ids = cell(1,1);

    if exist(fname, 'file')

        % load stepData
        load(fname, 'stepData');
        
        % load neural data
        load(fullfile(getenv('SSD'), 'paper2', 'modelling', 'neuralData', [session '_neuralData.mat']), ...
            'unit_ids', 'spkRates', 'timeStamps');
        t = timeStamps;  % for neural data
        nunits = length(unit_ids);
        
        % load predicted neural data
        predictedRates = nan(size(spkRates));
        for j = 1:length(unit_ids)
            file = fullfile(getenv('SSD'), 'paper2', 'modelling', 'glms', 'residual_glms', ...
                [session '_cell_' num2str(unit_ids(j)) '_glm.mat']);
            if exist(file, 'file')
                load(file, 'fitdata');
                predictedRates(j,:) = interp1(fitdata.t, fitdata.yhat, t);
            end
        end
        
        % compute responses
        for j = 1:4
            nsteps = height(stepData{j});
            resp = nan(nsteps, nunits, nx);
            pred = nan(nsteps, nunits, nx);

            for k = 1:nsteps
                epoch = stepData{j}.times(k,:);  % start and end times for step
                bins = t>=epoch(1) & t<epoch(2);
                if sum(bins)>1
                    tsub = t(bins);
                    if nunits==1
                        resp(k,:,:) = interp1(tsub, spkRates(bins), ...
                            linspace(tsub(1), tsub(end), nx), 'linear');
                        pred(k,:,:) = interp1(tsub, predictedRates(bins), ...
                            linspace(tsub(1), tsub(end), nx), 'linear');
                    else
                        resp(k,:,:) = interp2(tsub, (1:nunits)', spkRates(:,bins), ...
                            linspace(tsub(1), tsub(end), nx), (1:nunits)', 'linear');
                        pred(k,:,:) = interp2(tsub, (1:nunits)', predictedRates(:,bins), ...
                            linspace(tsub(1), tsub(end), nx), (1:nunits)', 'linear');
                    end
                end
            end
            responses{j} = permute(resp, [2 1 3]);  % unit X step X 'time' (after permutation)
            predicted{j} = permute(pred, [2 1 3]);  % unit X step X 'time' (after permutation)
        end
    end
    outs = {stepData, responses, predicted, unit_ids};  % pack outputs to allow parallelization
end

