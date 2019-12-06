function colors = colorme(numColors, varargin)

s.offset = .1;               % offset for phase of colors (0->1)
s.saturation = 1;            % saturation for all colors
s.value = 1;                 % value for all colors
s.showSamples = true;        % whether to make figure with sample of color scheme
s.bgColor = 'white';         % background color for sample figure

% reassign settings passed in varargin
if exist('varargin', 'var'); for i = 1:2:length(varargin); s.(varargin{i}) = varargin{i+1}; end; end

% get hue values
h = linspace(0,1,numColors+1);
h = h + s.offset/numColors;
h = rem(h, 2);

% make colors!
colors = nan(numColors,3);
for i = 1:numColors
    colors(i,:) = hsv2rgb([h(i) s.saturation s.value]);
end

% show color samples
if s.showSamples
    figure('name', sprintf('%i colors, %.2f offset', numColors, s.offset), ...
        'color', s.bgColor, 'menubar', 'none', 'position', [400 400 400 200]); hold on
    x = linspace(0,2*pi,100);
    phaseOffsets = linspace(0,2*pi,numColors+1);
    for i = 1:numColors
        plot(sin(x+phaseOffsets(i)), 'linewidth', 10, 'color', colors(i,:));
    end
    set(gca, 'visible', 'off')
end
