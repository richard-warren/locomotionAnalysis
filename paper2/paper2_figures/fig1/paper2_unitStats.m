% compute basic stats about number of units recorded
data = getUnitInfo('nucleiOnly', false);
data = data(~strcmp(data.nucleus, 'unregistered'), :);
%%

height(data)