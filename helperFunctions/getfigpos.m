function figPos = getfigpos()

figPos = get(gcf, 'Position');
figPosText = sprintf('[%i %i %i %i]', figPos);
clipboard('copy', figPosText)
disp(figPosText)
