function [k, p] = ktestnull(smp, varargin)

% compute kuiper statistic for difference between circular distribution of
% phase samples in smp, and circular uniform distribution // smphsould
% range between 0 and 2*pi // if 'bootstraps'>0, computes p via
% bootstrapping procedure: computes k stats for bootstraps random samples
% of length = length(smp), fits a gaussian to k distribution, and computes
% p as 1 - integral of gaussian to the left of the k stat (one sided test)

% todo: should be easy to modify such that two samples can be compared,
% rather than one sample being compared to uniform distribution! // should
% also check gaussian-ness of bootstrapped k distribution... is it safe to
% assume this is gaussian???



% settings
s.res = 40;
s.plot = false;
s.bootstraps = 0;


% inits
if exist('varargin', 'var'); for i = 1:2:length(varargin); s.(varargin{i}) = varargin{i+1}; end; end % reassign settings passed in varargin

% compute k stat for emperical distribution
[k, dpind, dmind, b, cdf, cdf_null] = getk(smp, s.res);


if s.bootstraps>0
    ks = nan(1, s.bootstraps);
    n = length(smp);
    for i = 1:s.bootstraps
        ks(i) = getk(rand(1,n)* 2*pi, s.res);
    end
    
    % fit gaussian
    p = 1 - normcdf(k, mean(ks), std(ks));
else
    p = nan;
end


if s.plot
    figure('color', 'white', 'menubar', 'none', ...
        'position', [400 100 381.00*((s.bootstraps>0)+1) 290.00])
    if s.bootstraps>0; subplot(1,2,1); end
    
    x = b(1:end-1) + (b(2)-b(1))/2;  % bin centers
    plot(x, [cdf; cdf_null]', '.-'); hold on
    plot([x(dpind), x(dpind)], [cdf(dpind), cdf_null(dpind)],'k-') % plotting vertical lines
    plot([x(dmind), x(dmind)], [cdf(dmind), cdf_null(dmind)],'k-') % plotting vertical lines
    title(sprintf('kuiper test: k=%.3f', k))
    set(gca, 'box', 'off')
    
    if s.bootstraps>0
        subplot(1,2,2); hold on
        
        xres = std(ks) / 100;
        xg = 0 : xres : max([ks+range(ks) k]);
        gaus = normpdf(xg, mean(ks), std(ks));

        histogram(ks, 20, 'Normalization', 'pdf')
        plot(xg, gaus, 'LineWidth', 3);
        set(gca, 'box', 'off', 'xlim', [min(xg(1),k) max(xg(end),k)], 'YColor', 'none')
%         set(gca, 'box', 'off', 'xlim', [xg(1) xg(end)], 'YColor', 'none')
        xlabel('bootstrap k distribution')
        ylims = ylim;
        plot([k k], ylims, 'color', 'black', 'LineWidth', 1)
        text(k, mean(ylims), sprintf('k=%.1e', k), 'HorizontalAlignment', 'right')
        title(sprintf('p=%.1e', p))
    end
end


% compute k statistic for a single sample
function [k, dpind, dmind, b, cdf, cdf_null] = getk(smp, res)
    
    % compute CDFs
    b = linspace(0, 2*pi, res+1);  % bin edges
    cdf = zeros(1, res);
    cdf_null = linspace(0, 1, res);

    for j = 1:res
        cdf(j) = sum(smp>=b(j) & smp<b(j+1));
    end
    cdf = cumsum(cdf);
    cdf = cdf / cdf(end);

    % get k statistic
    [dplus,  dpind] = max([0 cdf-cdf_null]);
    [dminus, dmind] = max([0 cdf_null-cdf]);
    k = dplus + dminus;
end
end




