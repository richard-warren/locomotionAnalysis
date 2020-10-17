function [dmat, t] = makeDesignMatrix(session, varargin)

% settings
s.timeDegrees = 0;  % add time polynomial of degree timeDegrees (constant term is excluded)

% inits
fprintf('making design matrix for %s... ', session)
if exist('varargin', 'var'); for i = 1:2:length(varargin); s.(varargin{i}) = varargin{i+1}; end; end  % parse name-value pairs
load(fullfile(getenv('SSD'), 'paper2', 'modelling', 'predictors', [session '_predictors.mat']), 'predictors');
settings = readtable(fullfile(getenv('GITDIR'), 'locomotionAnalysis', 'paper2', 'glm', 'predictorSettings.csv'));
settings = settings(logical(settings.include),:);
t = predictors.t{1};  % assumes first predictor is continuous
dt = t(2)-t(1);
rows = length(t);
dmat = table();


% add predictors one at time, apply settings defined in predictorSettings.csv
for i = 1:height(settings)
    name = settings.name{i};
    type = predictors{name, 'type'};
    data = predictors{name, 'data'}{1};
    
    if type=='event'
        if any(isnan([settings.kernel_start(i) settings.kernel_stop(i) settings.n_kernels(i)]))
            fprintf('WARNING! Settings missing for predictors %s\n', name)
            break
        end
        
        data = data(~isnan(data));
        binned = histcounts(data, t(1)-dt/2 : dt : t(end)+dt/2);
        kernels = makeCosBasis(settings.kernel_start(i), settings.kernel_stop(i), settings.n_kernels(i), 'dt', dt);
        bases = nan(length(t), settings.n_kernels(i));
        for j = 1:settings.n_kernels(i); bases(:,j) = conv(binned, kernels(j,:), 'same'); end
        addPredictor(bases, name)
        
    elseif type=='epoch'
        if any(isnan([settings.n_kernels(i) settings.max_epoch_duration(i)]))
            fprintf('WARNING! Settings missing for predictors %s\n', name)
            break
        end
        
        durations = diff(data,1,2);
        data = data(durations<settings.max_epoch_duration(i),:);
        bases = zeros(length(t), settings.n_kernels(i));
        
        for j = 1:size(data,1)
            pad = (data(j,2) - data(j,1)) *.25;  % pad window on sides so edges of first and last cosine bump fit
            bins = t>(data(j,1)-pad) & t<(data(j,2)+pad);
            if any(bins)
                basis = makeCosBasis(data(j,1), data(j,2), settings.n_kernels(i)+1, 't', t(bins));
                bases(bins,:) = bases(bins,:) + basis(1:end-1,:)';
            end
        end
        addPredictor(bases, name)
        
    elseif type=='continuous'
        addPredictor(data, name)
    end
end

% add time polynomial
bases = nan(rows, s.timeDegrees);
tNorm = t - mean(t);  % not sure if this affects anything...
for i = 1:s.timeDegrees
    bases(:,i) = tNorm.^i;
end
addPredictor(bases, 'time')




fprintf('all done!\n')



function addPredictor(data, name)
    % add predictor column to dmat (or initialize dmat if it is empty)
    
    if size(data,1)~=rows; data = data'; end  % orient correctly
    
    if isempty(dmat)
        dmat = table(data, 'VariableNames', {name});
    else
        dmat = [dmat table(data, 'VariableNames', {name})];
    end
end
end