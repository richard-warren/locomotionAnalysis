function [ radius, center ] = circfit( data, weights )
%CIRCFIT uses Weighted Linear Least Squares method to find the circle best
%fitting the data (x/y-measurements as columnvectors of the data-matrix).
%Circle fitting is done according to a method proposed by I.D. Coope in
%"Coope, I. D.,Circle Fitting by Linear and Nonlinear Least Squares,
%University of Canterbury, Mathematics Department, Report No. 69, 1992."
%and exteded to use Weighted Linear Least squares (idea and implementation
%by D.Buehler - buehler.dani@gmail.com).


%transform data according to Coopes method
A=[data;ones(1,size(data,2))]';
y=sum(data.^2,1)';

%Weighted Linear Least Squares
p=(A'*diag(weights)*A)\A'*diag(weights)*y;

%resubstitue to find desired parameters
center = .5 * p(1:end-1,:);
radius = sqrt(p(end,1)+center'*center);
end

