function thresh = getScoreThresh(session, metadata, varargin)
% determines the confidence (score) threshold for neural network tracking,
% which depends on whether deeplabcut or deepposekit was used // dlc uses
% binary output, and needs higher threshold, and dpk uses gaussian outputs,
% and needs lower threshold // this function reads metadata associated with
% each session to determine the appropriate threshold // metadata is the
% name of the metadata file

s.defaultThresh = .99;  % thresh if no metadata in session folder
s.deeplabcutThresh = .99;
s.deepposekitThresh = .8;
s.verbose = false;


% initializations
if exist('varargin', 'var'); for i = 1:2:length(varargin); s.(varargin{i}) = varargin{i+1}; end; end % reassign settings passed in varargin

metaFile = fullfile(getenv('OBSDATADIR'), 'sessions', session, metadata);

if exist(metaFile, 'file')
    load(metaFile, 'analysis')    
    
    if strcmp(analysis, 'deepposekit')
        thresh = s.deepposekitThresh;
        if s.verbose; fprintf('%s: using deepposekit confidence threshold: %.2f\n', session, thresh); end
    elseif strcmp(analysis, 'deeplabcut (old)')
        thresh = s.deeplabcutThresh;
        if s.verbose; fprintf('%s: using deeplabcut confidence threshold: %.2f\n', session, thresh); end
    end
else
    thresh = s.defaultThresh;
    if s.verbose; fprintf('%s: using default confidence threshold: %.2f\n', session, thresh); end
end