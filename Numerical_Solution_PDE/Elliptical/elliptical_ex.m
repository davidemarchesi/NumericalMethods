% The aim of this script is to test the function DirEllipt for the
% resolution of elliptical equations with the Finite Elements Method.
% It is taken a simplified 1D model of omogeneus Dirichlet equation.
% Will be examined also the convergence properties and the characteristics
% of the equations.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% -mu u'' + beta u' + cu = f(x), in [a,b],   %
% u(a) = u(b) = 0,                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
% In this part will be considered the equation with the following coeff. :
mu = 1 ;
beta = 2 ;
c = 3 ;
g =@(x) (x-x)-1 ;
a = 0 ;
b = 1 ;
h = 1/100 ;

% Call of the solving function :
[X,U] = DirEllipt(mu,beta,c,g,a,b,h) ;

% Definition of the exact solution :
C_2 = 1/3 * ((exp(1)-1)/(exp(4)-1)) ;
C_1 = 1/3 - C_2 ;
u_ex =@(x) (C_1.*exp(-x))+(C_2.*exp(3*x))-(1/3) ;

% Visualization of the numerical resolution vs the exact solution :
figure(1)
xlabel('x') ;
ylabel('u') ;
title('Numerical solution vs Exact solution');
hold on 
plot(X,U,"LineWidth",1.6,"Marker","diamond","Color",'red');
step = linspace(a,b,100) ;
plot(step, u_ex(step),"LineWidth",1.6,"Color",'blue') ;
legend('Num. Sol.','Ex. Sol.','Location','eastoutside');
hold off

% A better visualization is given by the loglog plane :
figure(2)
U_log = U(2:end,1);
X_log = X(2:end,1);
semilogx(X_log,U_log,"LineWidth",1.6,"Marker","diamond","Color",'red');
step = linspace(a+0.01,b,99) ;
hold on
semilogy(step, u_ex(step),"LineWidth",1.6,"Color",'blue') ;
xlabel('x') ;
ylabel('u') ;
title('Numerical solution vs Exact solution');
legend('Num. Sol.','Ex. Sol.','Location','eastoutside');
hold off
%%
% Now will be considered the equation with the following values :
mu = 1;
beta = 0;
c = 0;
z =@(x) -2.*exp(x).*cos(x) ;
a = 0 ;
b = 2*pi ;
h = 2*pi/60 ;

% The exact solution is :
uex =@(x) exp(x).*sin(x) ;

% Now we use the function :
[Xnew,Unew] = DirEllipt(mu,beta,c,z,a,b,h) ;

% Here the plot of the numerical and of the exact function :
figure(3)
xlabel('x');
ylabel('u');
title('Exact solution vs Numerical solution');
plot(Xnew,Unew,"LineWidth",1.6,"Color",'magenta','Marker','o');
hold on
q = linspace(a,b,60) ;
plot(q,uex(q),"LineWidth",1.6,"Color",'green');
hold off

% Now we want to do an evaluation of the error using the L^2 norm.
% In order to do this we have to previously calculate the integral of the
% linear function given by the solution of the equation.
% So we have two steps :
% 1- Interpolation of the function (linear interpolation)
w = linspace(a,b) ;
u_h =@(x) interp1(Xnew,Unew,x) ;
err =@(x) (uex(x)-u_h(x)).^2 ;
N = length(Xnew) ;
% 2- Integration of the function using the Simpson composite formula
int_err = simpson_int(a,b,N,err) ;

normL2 = sqrt(int_err) ;

% Doing a plot of the error :
figure(4)
xlabel('x') ;
ylabel('u') ;
title('Plot of exact fun., num. fun. and errors (^2)');
hold on
plot(w,u_h(w),"LineWidth",1.5,"LineStyle","--","Marker","*");
plot(w,uex(w),"LineWidth",1.5);
plot(w,err(w),"LineWidth",1.5,"LineStyle",":","Marker","o");
legend('Num. Sol.','Ex. Sol.','Error',Location='eastoutside');
hold off

% Now we want to extend this analysis to different values of 'h'. Comparing
% the different errors obtained :
norm_vec = [] ;
Nv = [20 ; 40; 80; 160; 320; 640] ;
figure(5)
for i = 1 : 1 : 6
    N = Nv(i);
    hnew = 2*pi/N ;
    [Xv,Uv]=DirEllipt(mu,beta,c,z,a,b,hnew) ;
    u_h =@(x) interp1(Xv,Uv,x) ;
    err =@(x) (uex(x)-u_h(x)).^2 ;
    int_err = simpson_int(a,b,N,err) ;
    normL2 = sqrt(int_err) ;
    norm_vec = [norm_vec ; normL2];
    loglog(w,err(w),'LineWidth',1.6);
    hold on
end
grid on
xlabel('x');
ylabel('Errors');
title('Plot of the errors with different step "h"');
legend('20','40','80','160','320','640',Location='eastoutside');
hold off
figure(6)
title('Norm L^2 in function of the step "h"');
xlabel('Number of subintervals');
ylabel('NormL2');
hold on
plot(Nv,norm_vec,'LineWidth',1.6);
hold off

return