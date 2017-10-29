function [x, y, scores] = nonMaximumSupress(scoresRaw, boxSize, thresh)

% temp
% scores = frameFiltered;
% boxSize = [36 36];
% thresh = .01;

% convert scores to 2D coordinates and sort
detected = find(scoresRaw(:) > 0);
[y, x] = ind2sub(size(scoresRaw), detected);
[scores, sortInds] = sort(scoresRaw(detected), 'descend');
y = y(sortInds);
x = x(sortInds);

% find box limits
x1 = y - .5 * boxSize(2);
x2 = y + .5 * boxSize(2);
y1 = x - .5 * boxSize(1);
y2 = x + .5 * boxSize(1);



% initializations
nPoints = length(detected);
area = prod(boxSize);
keepers = true(1, length(detected));

for i = 1:nPoints
    for j = (i+1):nPoints
        
        if keepers(j)
            
            xOverlap = max(0, min(x2(i), x2(j)) - max(x1(i), x1(j)));
            yOverlap = max(0, min(y2(i), y2(j)) - max(y1(i), y1(j)));

            overlap = xOverlap * yOverlap;
            union = 2*area - overlap;

            if( (overlap / union) > thresh)
                keepers(j) = false;
            end
        end
    end
end

scores = scores(keepers);
y = y(keepers);
x = x(keepers);
