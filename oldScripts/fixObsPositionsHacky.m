function obsPositionsFixed = fixObsPositionsHacky(obsPositions)

% !!! this code should be replaced with something that does a better job of determing the obs position
% relative to the nose without the deeplabcut analysis needing to be performed first

nosePosApprox = .325;
obsPositionsFixed = obsPositions - nosePosApprox;