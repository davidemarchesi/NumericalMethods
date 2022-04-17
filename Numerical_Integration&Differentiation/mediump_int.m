function I = mediump_int(a, b, N, fun)
% Integral Calculated with the medium point method
% Inputs:
% - a,b: Integration Extremes;
% - N: Number of subintervals;
% - fun: function;
% Output:
% - I: integral;

% Subintervals length
h = (b-a)/N;

% Subintervals extremes/medium points
x = (a+h/2):h:(b-h/2);

% Evaluation of the points
fx = fun(x);

% Approximated integral calculation
I = 0;
for j=1:N
    I=I+h.*fx(j);
end
end