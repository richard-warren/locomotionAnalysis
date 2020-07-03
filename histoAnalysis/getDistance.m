function d = getDistance(point1, point2)


substracted = point2 - point1; % direction is from point1 to point2
d = sqrt(sum(substracted.^2));

end