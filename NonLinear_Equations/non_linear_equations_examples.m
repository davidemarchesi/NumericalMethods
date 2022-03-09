%%
% This file is full of plot and exercises related with the
% numerical solutions for non-linear equations
%%
% Here there's an example of function on which the methods can be
% applied:
xx=linspace(-2,2,10^(6));
f=@(x) exp(x)- x.^2 - sin(x) -1;
plot(xx , f(xx), 'b');
hold on;
y=0*ones(size(xx));
plot(xx, y, 'r');
hold on;
w=0;
r=0;
plot(w,r,'k:*','LineWidth',1);
title('Function Example for the bisection method');
xlabel('Abscissas');
ylabel('Ordinates');
%%
% How to use the bisection method function:
f = @(x) exp(x)-x.^2-sin(x)-1;
a=-2;
b=2;
toll=10^-12;
nmax=10^6;
[xvect,xdif,fx,it]=bisection_method(a,b,nmax,toll,f);
%%
% How to use the Newton method function:
f=@(x) exp(x)-x.^2-sin(x)-1;
df=@(x) exp(x)-2.*x-cos(x);
tol=10^-12;
nmax=10^6;
x0= 1.7;
[ xvect, it ] = newton_method(x0 , nmax , tol , f , df );
