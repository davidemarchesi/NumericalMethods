function [t_h,u_h,iter_fp]=eulero_bwd(f,t0,t_max,y_0,h)
% This function is the 'Backward Eulero' numerical method to solve ODE , in
% particular the Cauchy Problems :
% y'=f(t,y)
% y(0)=y_0
%
% Using the Backward Eulero method :
% (Implicit Eulero) : u_(n+1)=u_n+h*f_(n+1).
% 
% [Obs.1 For each time instant, for u_(n+1) the equation (not linear a priori)
%  is solved with a fixed point method : u=u_n+delta_t*f(t_(n+1),u); 
%  This method uses the iterative function : phi(x)=u_n + h *f(t_(n+1),x)]
%
% [Obs.2 The convergence condition for the fixed point method (|phi(x)'|<1)
%  is : h*|\partial_x f(t_(n+1),x)| < 1 ;
%  This condition is granted when 'h' results to be sufficiently small]
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
% - iter_fp : Vector which contains the number of fixed point iterations
%             used to solve the non-linear equation at every time instant;

% Discretization
t_h = t0:h:t_max;
N=length(t_h);
u_h=zeros(1,N);

% Iterative cycle that calculates : u_(n+1)=u_n+h*f_(n+1) . 
% (Note that at everi time iteration, fixed point sub-iterations must be 
%  made to calculate u_(n+1));
 
u_h(1)=y_0;

% Parameters for the fixed point iteration using the function 'fixpt.m'
N_max=300;
tol=1e-5;
iter_fp=zeros(1,N);

for it = 2:N
    
    % Variables declaration for the sub-iterations
    u_old=u_h(it-1);
    t_fp=t_h(it);
    
    phi=@(u) u_old + h * f( t_fp, u );
    % Declaration of 'phi' (the variables of 'phi' are the ones in the 
    %                       parentesis, the others (all in exception of 'u')
    %                       are interpreted by Matlab as constants)
    % 
    % Subiterations
    [u_fp, it_fp] = fixpt(u_old, phi, N_max, tol);
    u_h(it)=u_fp(end);
    
    % Fixed point iterations (can be used to evaluate the convergence of the method)
    iter_fp(it)=it_fp;
    
end
end