function I = gauss_sq_I(a,b,fun,m)
% Squaring Gauss Formulas with interpolating polynomials of grade n=m-1
% Inputs:
% - a,b: Integral Extremes;
% - fun: Function;
% - m  : Number of nodes;
% Output:
% - I  : Integral;
% This is a simple script of the Gauss Squaring formula considering a
% maximum number of points (5), and as a consequence also a maximum grade
% of the interpolating polynomial (4)
if m>5
    error('Maximum number of points: m=5')
end
% Gauss nodes table in the Interval [-1,1]
xbar=[
    0 0 0 0 0; % m=1
    -1/sqrt(3) 1/sqrt(3) 0 0 0; % m=2
    -sqrt(3/5) 0 sqrt(3/5) 0 0; % m=3
    -sqrt(3/7+2/7*sqrt(6/5)) -sqrt(3/7-2/7*sqrt(6/5)) sqrt(3/7-2/7*sqrt(6/5)) sqrt(3/7+2/7*sqrt(6/5)) 0; % m=4
    -1/3*sqrt(5+2*sqrt(10/7)) -1/3*sqrt(5-2*sqrt(10/7)) 0 1/3*sqrt(5-2*sqrt(10/7)) 1/3*sqrt(5+2*sqrt(10/7)) % m=5
    ];
% Gauss wights table in the Interval [-1,1]
wbar=[
    2 0 0 0 0 % m=1
    1 1 0 0 0 % m=2
    5/9 8/9 5/9 0 0 % m=3
    (18-sqrt(30))/36 (18+sqrt(30))/36 (18+sqrt(30))/36 (18-sqrt(30))/36  0 % m=4
    (322-13*sqrt(70))/900 (322+13*sqrt(70))/900  128/225 (322+13*sqrt(70))/900  (322-13*sqrt(70))/900 % m=5
    ];
% Are selected the correspondent nodes to the grade n=m-1 related to the
% interpolating polynomial
%
% Transformation of the nodes in the interval [a,b]
xi=(a+b)/2+(b-a)/2*xbar(m,:);
%
% Trasformation of the weights in the interval [a,b]
wi=wbar(m,:)*(b-a)/2;
% Integral Calculation
I=fun(xi)*wi';
return