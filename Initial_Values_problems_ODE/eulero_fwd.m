function [t_h,u_h]=eulero_fwd(f, t_0, t_max, y_0, h)
% This function is the 'Forward Eulero' numerical method to solve ODE , in
% particular the Cauchy Problems :
% y'=f(t,y)
% y(0)=y_0
%
% Using the Forward Eulero method :
% (Explicit Eulero) : u_(n+1)=u_n+h*f_n
%
% Inputs :
% - f : Function which describes the Cauchy Problem . It must be formed by
%      2 variables : f=f(t,y);
% - t_0,t_max : Time extremes of the solution; 
% - y_0 : Starting datum of the Cauchy Problem;
% - h : length of the time discretization steps;
%
% Outputs :
% - t_h : Vector of the time at which the approximated solution is
%         calculated;
% - u_h : Discreet solution calculated at every instant 't_h';

% Discretization :
t_h = t_0:h:t_max ;
N = length(t_h) ;
u_h = zeros(1,N) ;

% Assignment of the starting value :
u_h(1) = y_0 ;

% Variable Declaration
u_old = 0;

% Iterative cycle 
for i = 2:N
    u_old = u_h(i-1);
    u_h(i) = u_old + h*f(t_h(i-1),u_old);
end
end