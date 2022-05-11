function[t_h, u_h] = eulero_fwd_vec (f, t_0, t_max, y_0, h)
% This function is the 'Forward Eulero' numerical method to solve ODE , in
% particular the Cauchy Problems :
% y'=f(t,y)
% y(0)=y_0
%
% Using the Forward Eulero method :
% (Explicit Eulero) : u_(n+1)=u_n+h*f_n
%
% With this method , the ordinary differential equation can be, more in general
% a vectorial one (So where y and f(t,y) are in the set R^d , where d = 1,2,...)
% (Note that if d = 1 the scalar case is obtained)
%
% Inputs :
% - f : Function which describes the Cauchy Problem . It must be formed by
%      2 variables : f=f(t,y) (with y a generic vector of length d);
% - t_0,t_max : Time extremes of the solution;
% - y_0 : Starting datum of the Cauchy Problem (column vector of length d);
% - h : length of the time discretization steps;
%
% Outputs :
% - t_h : Vector of the time at which the approximated solution is
%         calculated;
% - u_h : Discreet solution calculated at every instant 't_h' (Obs. u_h is rectangular
%         matrix of dimension 'd x N_istants');

% Control that y_0 is a column vector
[ny_0, my_0] = size(y_0);
if ny_0 < my_0
   warning('y_0 must be a column vector! (Converted into a column one)');
   y_0 = y_0';
end

% Discretization :
t_h = t_0:h:t_max ;
N = length(t_h) ;
d = length(y_0) ;
u_h = zeros(d,N) ;

% Assignment of the starting value :
u_h(:,1) = y_0 ;

% Variable Declaration
u_old = 0;

% Iterative cycle
for i = 2:N
    u_old = u_h(:, i-1);
    u_h(:, i) = u_old + h*f(t_h(i-1),u_old);
end
end
