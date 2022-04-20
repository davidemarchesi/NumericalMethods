function I = simpson_int(a, b, N, fun)
% Integral Calculated with the Trapezium method
% Inputs:
% - a,b: Integration Extremes;
% - N: Number of subintervals;
% - fun: function;
% Output:
% - I: integral;

% Subintervals length
h = (b-a)/N;

% Subintervals extremes/medium points
x = a : h : b; 

% Evaluation of the points
fx = fun(x);

% Approximated integral calculation
I = 0;
xm = 0;
for j=1:N
    xm = 0.5.* (x(j+1)+x(j));
    I = I + h./6 .*(fx(j) + fun(xm).*4 +fx(j+1));
end