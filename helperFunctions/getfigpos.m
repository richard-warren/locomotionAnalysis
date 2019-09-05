function figPos = getfigpos()

figPos = get(gcf, 'Position');
figPosText = sprintf('[%.2f %.2f %.2f %.2f]', figPos);
clipboard('copy', figPosText)
disp(figPosText)
