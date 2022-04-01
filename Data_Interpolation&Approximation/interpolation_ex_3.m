% This script has the aim to underline the difference between the
% selection of the nodes in the application of the interpolating method
% with the Lagrange Polynomials
%
% Will be considered the Runge's function approximation problem :
f =@(x) 1./(1+x.^2);
x = linspace(-5,5);
%
% Now we build the Lagrange Interpolating Polynomials :
figure();
hold on;
title('Runge real fun. vs Lagrange Int. pol.');
xlabel('x-values');
ylabel('y-values');
plot(x,f(x),'LineWidth',1.75,'Color','r');
% Here the plot of the real, and the interpolated function
for n = [5,10]
    xx = linspace(-5,5,n+1);
    yy = f(xx);
    kl = polyfit(xx,yy,n);
    pl = polyval(kl,x);
    if n == 5
        clr = 1;
    end
    if n == 10
        clr = 0.2;
    end
    plot(x,pl,'LineWidth',1.4,'Color',[0 clr clr]);
end
legend('Runge Real','Lag. pol. n=5','Lag. pol. n=10');
hold off;
% Here the plot of the error between the real and the interpolated fun.
figure();
title('Error value for the Interpolating Lag. Pol. with equispaced nodes');
xlabel('x-values');
ylabel('Abs. Error');
hold on;
for n = [5,10]
    xx = linspace(-5,5,n+1);
    yy = f(xx);
    kl = polyfit(xx,yy,n);
    pl = polyval(kl,x);
    err = abs(f(x)-pl);
    max_err = max(err);
    mid_err = mean(err);
    fprintf('The maximum error for the Lagrangian Int.%d with equispaced nodes is : %f \n',n,max_err);
    fprintf('The medium error for the Lagrangian Int.%d with equispaced nodes is : %f \n',n, mid_err);
    plot(x,err,'LineStyle','--','LineWidth',1.5);
end
legend('L5 error','L10 error');
hold off;
%%
% Now will be done the same analysis, but this time using the CGL nodes :
%
figure();
hold on;
title('Runge real fun. vs Lag. with CGL');
xlabel('x-values');
ylabel('y-values');
plot(x,f(x),'LineWidth',1.75,'Color','r');
% Firstly we have to calculate them in the interval I=[-1,1], and then  we
% have to transform them to fit in the interval [-5,5]
xcgl = [];
xc = 0;
xc_1 = 0;
for n = [5,10]
    for i = 1:n+1
        % This is the calculation in the interval [-1,1]
        xc_1 = -cos(pi*i/(n+1));
        % This is the transformation we do to bring them in the interval [-5,5]
        xc = ((-5+5)/2)+(((5+5)/2)*xc_1);
        xcgl = [xcgl, xc];
    end
    Ycgl = f(xcgl);
    kl_cgl = polyfit(xcgl,Ycgl,n);
    pl_cgl = polyval(kl_cgl,x);
    if n == 5
        clr = 1;
    end
    if n == 10
        clr = 0.2;
    end
    plot(x,pl_cgl,'LineWidth',1.4,'Color',[0 clr clr]);
end
legend('Runge Real','Lag.5 CGL','Lag.10 CGL');
hold off
% Lets examinate now the error :
figure();
title('Error value for the Interpolating Lag. Pol. with CGL nodes');
xlabel('x-values');
ylabel('Abs. error');
hold on
xcgl = [];
xc = 0;
xc_1 = 0;
for n = [5,10]
    for i = 1:n+1
        xc_1 = -cos(pi*i/(n+1));
        xc = ((-5+5)/2)+(((5+5)/2)*xc_1);
        xcgl = [xcgl, xc];
    end
    Ycgl = f(xcgl);
    kl_cgl = polyfit(xcgl,Ycgl,n);
    pl_cgl = polyval(kl_cgl,x);
    err_cgl = abs(f(x)-pl_cgl);
    max_err_cgl = max(err_cgl);
    mid_err_cgl = mean(err_cgl);
    fprintf('The maximum error for the Lagrangian Int.%d with CGL nodes is : %f \n',n,max_err_cgl);
    fprintf('The medium error for the Lagrangian Int.%d with CGL nodes is : %f \n',n, mid_err_cgl);
    plot(x,err_cgl,'LineWidth',1.5,'LineStyle','--');
end
legend('L5 CGL error','L10 CGL error');
hold off;