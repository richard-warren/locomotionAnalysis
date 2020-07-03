% INPUTS: 
% variable = column vector of numeric measurements
% groups = class matching each instance in variable (either 0 or 1)
%
function [p_valuet, hypothesis] = Ttestplot(variable,groups)
%
% This function Performs an independent sample t-test with unequal variance
% assumed --  much more reliable (Zimmerman and Zumbo, 1993).
% Utilises the Box-Cox transformation trick to obtain normally distributed
% data. Plots the data with options to modify the box plots i.e. color the
% boxplot and change the line colors of the boxplot. Produces relevent
% descriptive statistics in a table and a t-test report in APA format
%
% By Musalula Sinkala (2017)
% University of Zambia - School of Health Sciences: Biomedical Sciences
% Dept. smsinks@gmail.com
%
% REFERECES
%
% Welch, B. L., 1947.The generalization of "Student's" problem when several
% different population variances are involved. Biometrika,34(1?2),28?35.
%
% Zimmerman, D.W. and Zumbo, B.D., 1993. Canadian Journal of Experimental
% Psychology 47(3), p.523.
%
% Massey, F. J., 1951. The Kolmogorov-Smirnov Test for Goodness of Fit. Journal
% of the American Statistical Association. Vol. 46, No. 253, pp. 68?78.
%
% Box GEP, Cox DR., 1964. An Analysis of Transformations. J R Stat Soc Ser B.
% 26:211?52.
%
% Tukey JW., 1977. Exploratory Data Analysis, Reading. Mass Addison-Wesley Publ
% Co.
%
% Error checking and data cleaning
v = variable;
for ii = 1:length(v)
    if  v(ii) < 0 && isnumeric(v(ii)) && isreal(v(ii))
        error(['All values of the measurement parameter needs to be '...
            'a positive integers']);
    end
end

for jj = 1:length(groups)
    if groups(ii) < 0 || groups(ii) > 1 && isnumeric(groups(ii)) && ...
            round(groups(ii)) ~= groups(ii)
        error('All values of the group must be either 0 or 1')
        % the groups need to be either 0 (for the controls) or 1 (cases)
    end
end

% view the distribution of the data
figure
hold on
histogram(v(groups==0),'FaceColor','c');
histogram(v(groups==1),'FaceColor','r');
fprintf('\n')
title('\bf Raw Data Distribution','FontSize',15)
legend({'Group 0','Group 1'},'FontSize',11,'Location','NorthWest')
ylabel('\bf Frequency','FontSize',12)
hold off

% clean the data
Data = [v,groups];
Data(any(isnan(Data),2),:) = [];
Data(any(Data(:,1)>=(median(Data(:,1))+3*std(Data(:,1))),2),:)=[];
Data(any(Data(:,1)<=(median(Data(:,1))-3*std(Data(:,1))),2),:)=[];
v = Data(:,1);
groups = Data(:,2);

[Shypothesis0,p_ktest0] = kstest(v(groups==0)); % Kolmogorov-Smirnov Test
[Shypothesis1,p_ktest1] = kstest(v(groups==1)); 
% group0 = p_ktest0
% group1 = p_ktest1

if Shypothesis1 == 0 || Shypothesis0 == 0
    fprintf('The data are normally distributed')
    % Plot the data
    figure
    hold on
    histogram(v(groups==0),'FaceColor','c');
    histogram(v(groups==1),'FaceColor','r');
    title('\bf Cleaned Data Distribution','FontSize',15)
    ylabel('\bf Frequency','FontSize',14)
    legend({'Group 0','Group 1'},'FontSize',11,'Location','NorthWest')
    hold off
    fprintf('\n');
    [hypothesis,p_valuet,ci,tstats] = ttest2(v(groups==0),v(groups==1)...
        ,'Vartype','unequal','Alpha',0.05);
    
else
    fprintf('The data are NOT normally distributed \n\n')
    % Requires the financial toolbox, Modify the hypothesis to invalid in
    % line and to either 0 or 1 in line 31 to switch OFF this block
    [transdat,lambda] = boxcox((v+1)); %(Box and Cox, 1964)
    quantile3 = quantile(transdat,3); InterQR = iqr(transdat);
    lowerfence = quantile3(1,1)-(InterQR*1.5);
    upperfence = quantile3(1,3)+(InterQR*1.5); %(Tukey, 1977)
    Data2 = [transdat,groups];
    Data2(any(Data2(:,1)<(lowerfence),2),:)=[];
    Data2(any(Data2(:,1)>(upperfence),2),:)=[];
    transdat = Data2(:,1);
    groups = Data2(:,2);
    
    figure
    hold on
    hist = histogram(transdat(groups==0),'FaceColor','c');
    histogram(transdat(groups==1),'FaceColor','r');
    title('\bf Box-Cox Transformed Data Distribution','FontSize',15)
    ylabel('\bf Frequency','FontSize',14)
    text(mean(transdat),max(hist.Values), sprintf('Lambda = %.5f',...
        lambda),'FontSize',15)
    legend({'Group 0','Group 1'},'FontSize',11,'Location','NorthWest')
    hold off
    [hypothesis,p_valuet,ci,tstats] = ttest2(transdat(groups==0),transdat...
        (groups==1),'Vartype','unequal','Alpha',0.05);
    
end

% T-test Stastics
T_test_stats = tstats
tst = getfield(T_test_stats,'tstat');
df = getfield(T_test_stats,'df');
upperci = ci(2,1);
lowerci = ci(1,1);
civ = (upperci+lowerci)/2;

% plot figure
figure
hold on
if Shypothesis1 == 0 || Shypothesis0 == 0 
    % this can be changed to 0 || 1 if No Financial Toolbox
    h = boxplot(v,groups);
else
    h = boxplot(transdat,groups);
end
set(h,'LineWidth',2.5);
set(gca, 'XTick',1:2, 'XTickLabel', {'Group 0','Group 1'},'FontSize',14)
title('\bf Mean Difference Across Groups','Fontsize',18)

% Change the boxplot color
a = get(get(gca,'children'),'children');   % Get the handles of all the objects
t = get(a,'tag');   % List the names of all the objects
box1 = a(5:6);   % try single for vector numbers to select
% different aspects of the box plot
set(box1, 'Color', 'b');% Set the color of the first box to green

% change the box plot fill colour
h = findobj(gca,'Tag','Box');
for j=1:length(h)
    patch(get(h(j),'XData'),get(h(j),'YData'),'y','FaceAlpha',.1);
end

% show the p-value in boxplot
if p_valuet > 0.001
    text(1.3,max(transdat), sprintf('p = %.4f',p_valuet),'FontSize',15);
else
    text(1.3,max(transdat), sprintf('p < 0.001'),'FontSize',15);
end
ylabel('\bf Levels','FontSize',12)
hold off

% Create a table of group statistics
Groups = {'Group 0';'Group 1'};
N = [length(v(groups==0));length(v(groups==1))];
Mean = [mean(v(groups==0));mean(v(groups==1))];
Median = [median(v(groups==0));median(v(groups==1))];
Range = [range(v(groups==0));range(v(groups==1))];
Std = [std(v(groups==0));std(v(groups==1))];
Min = [min(v(groups==0));min(v(groups==1))];
Max = [max(v(groups==0));max(v(groups==1))];
Descriptive_Stats = table(N,Mean,Median,Range,Std,Min,Max,'RowNames',Groups)
fprintf('\n')

% Report the t-test results APA format: Final report requires editing
fprintf('T-TEST REPORT\n')

fprintf('An Independent sample t-test was performed to determine')
fprintf(' if there were \ndifferences in mean XX levels between ')
fprintf('Group 0 and Group 1.\n')
if p_valuet < 0.05 && mean(v(groups==0)) >= mean(v(groups==1))
    fprintf('XX levels were higher in Group 0 ')
    fprintf(['(%2.3f ± %2.3f) than \nGroup 1 (%2.3f ± %2.3f). A '...
        'statistically significant difference \n'],mean(v(groups==0)),...
        std(v(groups==0)),mean(v(groups==1)),std(v(groups==1)))
    fprintf(['of %2.3f (95%% CI, %2.3f to %2.3f), t(%2.3f) = %2.3f',...
        ', p = %2.4f.\n'],civ,lowerci,upperci,tst,df,p_valuet);
elseif p_valuet < 0.05 && mean(v(groups==0)) <= mean(v(groups==1))
    fprintf('XX levels were higher in Group 1 \n')
    fprintf(['(%2.3f ± %2.3f) than \nGroup 0 (%2.3f ± %2.3f). A '...
        'statistically  significant difference'],mean(v(groups==1)),...
        std(v(groups==1)),mean(v(groups==0)),std(v(groups==0)))
    fprintf(['of %2.3f (95%% CI, %2.3f to %2.3f),\n t(%2.3f) = %2.3f',...
        ', p = %2.4f.\n'],civ,lowerci,upperci,tst,df,p_valuet);
else
    fprintf('XX levels were not different in Group 0 ')
    fprintf(['(%2.3f ± %2.3f) and \nGroup 1 (%2.3f ± %2.3f). A '...
        'non-statistically significant difference \n'],mean(v(groups==0)),...
        std(v(groups==0)),mean(v(groups==1)),std(v(groups==1)))
    fprintf(['of %2.3f (95%% CI, %2.3f to %2.3f), t(%2.3f) = %2.3f',...
        ', p = %2.4f.\n'],civ,lowerci,upperci,tst,df,p_valuet);
end
end

