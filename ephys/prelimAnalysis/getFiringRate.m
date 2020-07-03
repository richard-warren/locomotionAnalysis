function [spikeRate, times] = getFiringRate(spkTimes, varargin)

% settings
s.kernel = 'doubleExp';  % 'gauss', or 'doubleExp'
s.kernelSig = .02;       % (s) if a gaussian kernel is used
s.kernelRise = .005;     % (s) rise for double exponential kernel
s.kernelFall = .02;      % (s) fall for double exponential kernel

s.fs = 1000;             % (smps/s) sampling frequency of instantaneous firing rate
s.tLims = [];           % (s) [min max] times for x axis


% initializations
if exist('varargin', 'var'); for i = 1:2:length(varargin); s.(varargin{i}) = varargin{i+1}; end; end  % parse name-value pairs


% create kernel
if strcmp(s.kernel, 'gauss')
    kernelX = 0 : 1/s.fs : s.kernelFall*10;
    kernelX = [-fliplr(kernelX(2:end)) kernelX];  % ensures x axis is centered at 0
    kernel = arrayfun(@(x) (1/(s.kernelSig*sqrt(2*pi))) * exp(-.5*(x/s.kernelSig)^2), kernelX);
    kernel = kernel/sum(kernel);
    
elseif strcmp(s.kernel, 'doubleExp')
    kernelX = 0 : 1/s.fs : s.kernelFall*20;
    kernel = arrayfun(@(x) exp(-x/(s.kernelFall))-exp(-x/(s.kernelRise)), kernelX);
    kernel = [zeros(1,length(kernel)-1) kernel];  % -1 because kernel should have odd number of entries (to be centered at 0)
    kernel = kernel/sum(kernel);
end


% make time axis
if isempty(s.tLims)
    s.tLims = [min(spkTimes)-length(kernel)/s.fs, max(spkTimes)+length(kernel)/s.fs];
end
times = s.tLims(1) : (1/s.fs) : s.tLims(2);


% convert spike times to binary vector
binEdges = [times times(end)+(1/s.fs)] - (1/s.fs)*.5;
spikesVec = histcounts(spkTimes, binEdges);
spikesVec = spikesVec * s.fs;


% convolve
spikeRate = conv(spikesVec, kernel, 'same');


