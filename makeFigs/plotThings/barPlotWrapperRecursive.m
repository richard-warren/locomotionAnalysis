% function dvMatrix = barPlotWrapperRecursive(data, dv, varsToGet, varsGotten, varLevels, varLevelNames, varsToAvg)


data = flat;
dv = 'isTrialSuccess';
vars = {'isLightOn', 'condition'};
varLevels = {[0,1], {'saline', 'muscimol'}};
varLevelNames = {{'light off', 'light on'}, {'mus', 'sal'}};
varsToAvg = {'trial', 'session', 'mouse'};