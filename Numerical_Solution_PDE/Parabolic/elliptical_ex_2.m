% The aim of this script is to solve and analyze different type of partial
% differential equations with the functions implemented in this repository.
% In particular will be analized an elliptical function with Nuemann and
% Dirichlet conditions, and then a more complex function : the parabolic
% one, 1 dimensional, in function of a space variable 'x' and of a time
% variable 't'.
%
% (Surely these methods could be extended to multi-dimensional problems, and
% also to problems with different boundary conditions [Dirichlet, Neumann,
% Robin, Mixed,...])
%
%%
% In this first part will be analyzed the elliptical equation :
%
%  -u''+u' = cos(x) ,      x in (0,2pi);
%   u(0) = 0 ;
%   u'(2pi) = 0 ;
%
% The discussion of the resolution of this type of equation was already
% made in previous script, here will be seen only the resolution with the
% function 'Dir_Neuelliptic.m'.

% Declaration of the values in the function :
% 1 - Coefficients
mu = 1;
beta = 1;
c = 0;
% 2 - Function handle
fun =@(x) cos(x);
% 3 - Extremes
a = 0;
b = 2*pi;
% 4 - Neumann right extreme condition
gN = 0;
% 5 - Step Length (will be seen the solutions at different step lengths 'h')
h = [2*pi/5 ; 2*pi/10; 2*pi/20; 2*pi/30; 2*pi/60];

figure(1)
title('Numerical Solution with different step "h"');
xlabel('X');
ylabel('Numerical Solution');
hold on
for i = 1 : length(h)
    [X,U] = Dir_Neuelliptic(mu,beta,c,fun,a,b,gN,h(i));
    plot(X,U,"LineWidth",1.5);
    hold on
end
legend(Location="eastoutside");
hold off

%%
% Now will be seen the resolution with the implemented function 'Dir_parabolic'
% of a parabolic function. Here the function that will be considered :
%
%   u_t - mu u_xx = 0 ,      x in [0,1], t in [0,0.2];
%   u(0,t) = u(1,t) = 0 ,    t in [0,0.2];
%   
%            | -x            x in [0,1/5);
%   u(x,0) = |  x-2/5        x in [1/5,7/10);
%            |  1-x          x in [7/10,1];
%
% Obs : for different values of theta I have different resoluters :
%       - theta = 0 , explicit Eulero, which needs the stability condition
%         to be satisfied ;
%       - theta = 0.5 , Crank-Nicolson method (implicit with quadratic
%         convergence);
%       - theta = 1 , implicit Eulero (implicit, with linear convergence);
%
% Now will be observed different type of solutions :

% Declaration of the equation values :
mu = 0.5;
beta = 0;
c = 0;
fun =@(x,t) 0.*x+0.*t ;
g =@(x) -x .* (x < 1/5) +...
       + x-2/5 .* (x >= 1/5) .* (x < 7/10) +...
       + (1-x) .* (x >= 7/10) ;
a = 0;
b = 1;
t0 = 0;
t_max = 0.2;

% Now let's impose the following values and solve the equation with all the
% three different methods :
theta1 = 0;
theta2 = 0.5;
theta3 = 1;
h = 0.1;
tau = 0.01;

% Implicit Eulero solution :
[X3,U3,T3,dt_max3] = Dir_parabolic(mu,beta,c,fun,g,a,b,t0,t_max,h,tau,theta3);
figure(1)
title('Implicit Eulero Solution');
xlabel('X');
ylabel('Numerical Solution');
hold on
plot(X3,U3,LineWidth=1.4);
hold off

% Crank-Nicolson solution :
[X2,U2,T2,dt_max2] = Dir_parabolic(mu,beta,c,fun,g,a,b,t0,t_max,h,tau,theta2);
figure(2)
title('Crank-Nicolson Solution');
xlabel('X');
ylabel('Numerical Solution');
hold on
plot(X2,U2,LineWidth=1.4);
hold off

% Explicit Eulero solution :
[X1,U1,T1,dt_max1] = Dir_parabolic(mu,beta,c,fun,g,a,b,t0,t_max,h,tau,theta1);
figure(3)
title('Explicit Eulero Solution');
xlabel('X');
ylabel('Numerical Solution');
hold on
plot(X1,U1,LineWidth=1.4);
hold off

% Can be easily seen that using the Implicit eulero method we have a
% problem with the stability, while with the other methods, being implicit,
% the stability is ensured. 
% We can easily see how the solution with the explicit method changes and
% reaches stability with the step 'tau' doing some tests :

tauvec = [0.01; 0.005; 0.0025];

for i = 1 : length(tauvec)
    [Xe,Ue,Te,dt_maxe] = Dir_parabolic(mu,beta,c,fun,g,a,b,t0,t_max,h,tauvec(i),theta1);
    figure()
    title('Explicit Eulero Solution with different "delta t"');
    xlabel('X');
    ylabel('Numerical Solution');
    hold on
    plot(Xe,Ue,LineWidth=1.4);
    hold off
end
% (This is an important analysis as the explicit method has a minor 
%  computational cost, but we have to ensure for it the proper stability)
%%
% In the following lines will be visualized the resolution of what seen
% above with a 3D graphic
mu = 0.5;
beta = 0;
c = 0;
fun =@(x,t) 0.*x+0.*t ;
g =@(x) -x .* (x < 1/5) +...
       + x-2/5 .* (x >= 1/5) .* (x < 7/10) +...
       + (1-x) .* (x >= 7/10) ;
a = 0;
b = 1;
t0 = 0;
t_max = 0.2;
theta = 0;
h = 0.1;
tau = 0.001;
[X,U,T,dt_max] = Dir_parabolic(mu,beta,c,fun,g,a,b,t0,t_max,h,tau,theta);
figure("Name","Animated Solution");
title('3D visualization of the solution');
[XX,TT] = meshgrid(X,T);
for i = 2:length(T)
    surf(XX(1:i,:),TT(1:i,:),U(:,1:i)')
    ylim([0;0.2]);
    xlim([0;1]);
    xlabel('Length');
    ylabel('Time');
    zlabel('Numerical Solution');
    drawnow()
    pause(0.05);
end