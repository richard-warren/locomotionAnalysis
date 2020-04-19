function whiskerTrimSchematic()

% settings
centers = 0:3:3*5;
colors = repmat([51 204 255]/255, 6, 1) .* fliplr(linspace(0,1,6))';

faceWid = 1;
faceHgt = 2.2;
padWid = .4;
padLoc = [.35 1.4];
noseWid = .2;
noseHgt = .15;
wiskAngles = [.2 1] * (pi/2);
wiskLengths = [.8 1.5];
eyePos = [.6 .4];
eyeRad = .12;
furColor = hsv2rgb([.08 .7 .6]);
whiskerThickness = .5;

padX = [2 3 4 5, 1 2 3 4 5, 1 2 3 4 5, 1 2 3 4 5, 1 2 3 4 5];
padY = [1 1 1 1, 2 2 2 2 2, 3 3 3 3 3, 4 4 4 4 4, 5 5 5 5 5];



% conditions
whiskers = {{'\alpha', '\beta', '\gamma', '\delta', 'A1', 'B1', 'C1', 'D1', 'E1', 'A2', 'B2', 'C2', 'D2', 'E2', 'A3', 'B3', 'C3', 'D3', 'E3', 'A4', 'B4', 'C4', 'D4', 'E4'}, ...
    {'\alpha', '\beta', '\gamma', '\delta', 'A1', 'B1', 'C1', 'D1', 'E1', 'A2', 'B2', 'C2', 'D2', 'E2', 'A3', 'B3', 'C3', 'D3', 'E3', 'A4', 'B4', 'C4', 'D4', 'E4'}, ...
    {'\gamma', '\delta', 'A1', 'B1', 'C1', 'D1', 'E1', 'A2', 'B2', 'C2', 'D2', 'E2', 'A3', 'B3', 'C3', 'D3', 'E3'}, ...
    {'\gamma', '\delta', 'D1', 'E1', 'D2', 'E2', 'D3', 'E3'}, ...
    {'\delta', 'E1', 'E2', 'E3'}, ...
    {'\delta'}};
bilateral = [1 0 0 0 0 0];
whiskerColors = jet(length(whiskers{1}))*.9;



% plot mice
figure('menubar', 'none', 'color', 'white', 'position', [98 795 1061 242]); hold on
for i = 1:length(centers)
    plotWhiskers(whiskers{i}, centers(i), bilateral(i), colors(i,:))
end
daspect([1 1 1])
set(gca, 'Visible', 'off')
saveas(gcf, fullfile(getenv('OBSDATADIR'), 'papers', 'hurdles_paper1', 'figures', 'matlabFigs', 'whiskerTrimSchematic'), 'svg');

% add schematic
figure('menubar', 'none', 'color', 'white', 'position', [886 955 273 235]); hold on
scatter(padX, padY, 300, whiskerColors, 'filled')
text(padX, padY, whiskers{1}, 'Color', 'white', ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'FontWeight', 'bold')
set(gca, 'Visible', 'off')
saveas(gcf, fullfile(getenv('OBSDATADIR'), 'papers', 'hurdles_paper1', 'figures', 'matlabFigs', 'whiskerTrimPadSchematic'), 'svg');



function plotWhiskers(whiskers, center, bilateral, conditionColor)

    % pad layout
    pad = {'\alpha', '\beta', '\gamma', '\delta', ...
        'A1', 'B1', 'C1', 'D1', 'E1', ...
        'A2', 'B2', 'C2', 'D2', 'E2', ...
        'A3', 'B3', 'C3', 'D3', 'E3', ...
        'A4', 'B4', 'C4', 'D4', 'E4'};

    % initializations
    padX = padX - mean(padX);
    padY = padY - mean(padY);
    bins = ismember(pad, whiskers);
    whiskerColorsMasked = whiskerColors;
    if any(~bins)
        whiskerColorsMasked(~bins,:) = repmat([.15 .15 .15], sum(~bins), 1);
    end

    % face outline
    t = 0:0.01:pi;
    x = faceWid*cos(t);
    y = faceHgt*sin(t);
    patch(x+center, y, furColor, 'EdgeColor', 'none', 'LineWidth', 3, 'facealpha', .75)
    plot(x+center, y, 'color', 'black', 'LineWidth', 2)
    plot([-1 1]*faceWid+center, [0 0], 'color', conditionColor, 'LineWidth', 5)

    % nose, eyes
    t = 0:0.01:2*pi;
    x = noseWid*cos(t);
    y = noseHgt*sin(t) + faceHgt-.5*noseHgt;
    patch(x+center, y, 'black')
    rectangle('position', [eyePos(1)-eyeRad+center, eyePos(2)-eyeRad, 2*eyeRad, 2*eyeRad], 'curvature', [1 1], 'facecolor', 'black', 'edgecolor', 'none');
    rectangle('position', [-eyePos(1)-eyeRad+center, eyePos(2)-eyeRad, 2*eyeRad, 2*eyeRad], 'curvature', [1 1], 'facecolor', 'black', 'edgecolor', 'none');

    % whisker pad
    scaling = (padWid/range(padX));
    x = padX*scaling + padLoc(1);
    y = padY*scaling + padLoc(2);
    scatter(x+center, y, 5, whiskerColorsMasked, 'filled');
    if bilateral
        scatter(-x+center, y, 5, whiskerColorsMasked, 'filled');
    else
        scatter(-x+center, y, 5, [.15 .15 .15], 'filled');
    end

    % whiskers
    xFrac = (padX-min(padX)) / range(padX);
    yFrac = (padY-min(padY)) / range(padY);
    for w = find(bins)
        a = interp1([-1 1], wiskAngles, -xFrac(w) + .1*yFrac(w));
        l = interp1([-1 1], wiskLengths, xFrac(w)-yFrac(w));
        plot([x(w) x(w)+cos(a)*l]+center, [y(w) y(w)+sin(a)*l], ...
            'Color', [.15 .15 .15], 'linewidth', whiskerThickness)
        if bilateral
            plot(-[x(w) x(w)+cos(a)*l]+center, [y(w) y(w)+sin(a)*l], ...
                'Color', [.15 .15 .15], 'linewidth', whiskerThickness)
        end
    end 
end
end