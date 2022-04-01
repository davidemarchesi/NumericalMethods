% The aim of this script is to compare the main interpolation models, which
% are: Lagrange's polynomials, Linear Composite, Natural Cubic Spline,
% Minimum Squared polynomial of grade 1-2-4.
%
% In order to do this we will assume to have done an experimental research
% about the link formula between the stress 's' and the deformation 'd' of
% a material, with the following values collected :
s = [0.00, 0.06, 0.14, 0.25, 0.31, 0.47, 0.60, 0.70]; %stress
d = [0.00, 0.08, 0.14, 0.20, 0.23, 0.25, 0.28, 0.29]; %deformation
%%
% Lagrange polynomial interpolation :
figure();
hold on;
plot(s,d,'*','LineWidth',1.5,'Color','r');
title('Real values vs Lagrange Polynomials');
xlabel('x-values');
ylabel('y-values');
n = length(s);
xx = linspace(0,0.8);
K = polyfit(s,d,n-1);
P = polyval(K,xx);
plot(xx,P,'LineWidth',1.5,'Color','b');
legend('Experimental values','Lagrange polynomial');
hold off;
%%
% Linear composite Interpolation :
figure();
hold on;
plot(s,d,'*','LineWidth',1.5,'Color','r');
title('Real values vs Linear Composite');
xlabel('x-values');
ylabel('y-values');
n = length(s);
xx = linspace(0,0.8);
% (the function 'interp1' takes three vectors, the first two are the ones 
% individuate the interpolating nodes; the last is the one which tells to 
% 'interp1' where to calculate the the function)
Comp_L = interp1(s,d,xx);
plot(xx,Comp_L,'LineWidth',1.5,'Color','c');
legend('Experimental values','Linear Composite');
hold off;
%%
% Natural Cubic Spline Interpolation :
figure();
hold on;
plot(s,d,'*','LineWidth',1.5,'Color','r');
title('Real values vs Cubic Spline');
xlabel('x-values');
ylabel('y-values');
n = length(s);
xx = linspace(0,0.8);
% Is necessary the function 'nat_spline'
Spline = nat_spline(s,d,xx);
plot(xx,Spline,'LineWidth',1.5,'Color','m');
legend('Experimental values','Cubic Spline');
hold off;
%%
% Minimum squared interpolation :
figure();
hold on;
plot(s,d,'*','LineWidth',1.5,'Color','r');
title('Real values vs Minimum squared');
xlabel('x-values');
ylabel('y-values');
n = length(s);
xx = linspace(0,0.8);
% obs. that the functions polyfit and polyval use the Lagrange polynomials
% only if the grade 'nn' inserted is equal to N-1, where N is the number of
% nodes. If the value inserted is < of N, then the function will use the
% minimum squared method.
for nn=[1,2,4]
    KK = polyfit(s,d,nn);
    MS = polyval(KK,xx);
    colr = 0.2*nn;
    plot(xx,MS,'LineWidth',1.5,'Color',[0+colr/2 1-colr 0+colr]);
end
legend('Experimental values','MS1','MS2','MS4');
hold off;
%%
% Now we want to compare all the methods with the experimental data in only
% one single plot :
figure();
hold on
title('Experimental data vs Interpolated functions');
xlabel('x-values');
ylabel('y-values');
xxx = linspace(0,0.7);
P2 = polyval(K,xxx);
Comp_L2 = interp1(s,d,xxx);
Spline2 = nat_spline(s,d,xxx);
plot(xxx,P2,'Color','g');
plot(xxx,Comp_L2,'Color','b');
plot(xxx,Spline2,'Color','c');
hold on;
for nnn=[1,2,4]
    KKK = polyfit(s,d,nnn);
    MSnew = polyval(KKK,xxx);
    plot(xxx,MSnew);
end
hold on
plot(s,d,'*','LineWidth',1.5,'Color','r');
legend('Lagrange P.','Lin. Comp.','Cub. Spline','MS1','MS2','MS4','Exp. Values');