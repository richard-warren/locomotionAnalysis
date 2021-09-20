function r2 = getModelR2(model, fitdata)
% model and fitdata produced by function like fitSequentialGlms.m

dt = fitdata.t(2) - fitdata.t(1);
yhat = exp(model.fit_preval) / dt;
y = fitdata.yRate;
r2 = corr(y', yhat)^2;
