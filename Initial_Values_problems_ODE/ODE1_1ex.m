% ODE exercises : this script has the aim to see the characteristics of
% some basic algorithms implemented to solve Ordinary Differential
% Equations, in particular Cauchy Problems. Will be seen the 'Forward Eulero' 
% algortihm and the 'Backward Eulero' algorithm.
%%
% Let's consider the following Cauchy Problem :
fun =@(t,y) (1./(1+t.^2))-2.*y.^2 ;
y_0 = 0;
% In the interval [0,10]
t_0 = 0;
t_max = 10;
t_vec = linspace(t_0,t_max) ;
% The analythic exact solution is :
y_ex =@(t) t./(1+t.^2) ;
%
%%
% Now we want to solve the function with the 'Forward Eulero' algorithm 
% using the implemented function 'eulero_fwd.m'.
%
% Time steps declaration :
t_h_fwd = [];
% Approximated outputs declaration :
u_h_fwd = [];
% Will be watched the solution for three different values of the step 'h' used
h = [0.2;0.1;0.05];
% Using the algorithm we obtain :
figure()
xlabel('time');
ylabel('function');
title('Calculated functions with Eulero Forward');
for k = 1:3
    step = h(k);
    [t_h_fwd,u_h_fwd] = eulero_fwd(fun,t_0,t_max,y_0,step) ;
    hold on
    plot(t_h_fwd,u_h_fwd,'LineWidth',1.4);
end
% To have a better visualization we also plot the analythic exact solution
hold on
plot(t_vec,y_ex(t_vec),'LineWidth',1.4,'LineStyle','-.');
legend('h=0.2','h=0.1','h=0.05','Exact solution');
hold off
%%
% Now let's visualize for the same problem a plot of the errors made by the
% Numerical method with different step lengths
%
% Consider the following steps :
h_err = [0.2,0.1,0.05,0.025,0.0125,0.00675];
%
% Declaration of the variables :
t_he_b = [];
u_he_b = [];
err_vec_b = [];
% Do the same iteration as seen above but this time calculating the errors
% step by step :
figure()
for d = h_err
    [t_he_b,u_he_b] = eulero_fwd(fun,t_0,t_max,y_0,d);
    % Calculation of the error at every point 't_he'
    err_vec_b = abs(y_ex(t_he_b) - u_he_b);
    % Calculation of the maximum error with the legth h_i
    err_max = max(err_vec_b);
    fprintf('\n The maximum error for the length %d (FW eu) , is : %f \n', d, err_max);
    % Logarithmic plot of the error made by the method :
    loglog(t_he_b,err_vec_b,'LineWidth',1.4);
    hold on;
end
grid on
xlabel('time');
ylabel('Absolute Error');
title('Comparison of the errors taking different time discretizations (FW eu)');
legend('0.2','0.1','0.05','0.025','0.0125','0.00675',Location='eastoutside');
hold off
%%
% Now in the following part of the script will be made the same
% considerations using the 'Backward Eulero' algorithm, and the
% corresponding implemented function 'eulero_bwd.m'
%
t_h_bwd = [];
u_h_bwd = [];
h = [0.2;0.1;0.05];
figure()
xlabel('time');
ylabel('function');
title('Calculated functions with Eulero Backward');
for k = 1:3
    step = h(k);
    [t_h_bwd,u_h_bwd,iter_fp] = eulero_bwd(fun,t_0,t_max,y_0,step) ;
    hold on
    plot(t_h_bwd,u_h_bwd,'LineWidth',1.4);
end
hold on
plot(t_vec,y_ex(t_vec),'LineWidth',1.4,'LineStyle','-.');
legend('h=0.2','h=0.1','h=0.05','Exact solution');
hold off
%%
% Error considerations for the 'Backward Eulero'
h_err = [0.2,0.1,0.05,0.025,0.0125,0.00675];
t_he_b = [];
u_he_b = [];
err_vec_b = [];
figure()
for d = h_err
    [t_he_b,u_he_b] = eulero_bwd(fun,t_0,t_max,y_0,d);
    err_vec_b = abs(y_ex(t_he_b) - u_he_b);
    err_max = max(err_vec_b);
    fprintf('\n The maximum error for the length %d (BW eu) , is : %f \n', d, err_max);
    loglog(t_he_b,err_vec_b,'LineWidth',1.4);
    hold on;
end
grid on
xlabel('time');
ylabel('Absolute Error');
title('Comparison of the errors taking different time discretizations (BW eu)');
legend('0.2','0.1','0.05','0.025','0.0125','0.00675',Location='eastoutside');
hold off
%%
% Now a comparison between the fwd and bwd solutions for a step h=1 :
h = 1;
[t_h_b2,u_h_b2,it_b2] = eulero_bwd(fun,t_0,t_max,y_0,h);
[t_h_f2,u_h_f2] = eulero_fwd(fun,t_0,t_max,y_0,h);
figure()
title('Two methods comparison for h=1');
xlabel('time');
ylabel('function values');
hold on
plot(t_h_f2,u_h_f2,'LineWidth',1.4,'LineStyle','--');
plot(t_h_b2,u_h_b2,'LineWidth',1.4,'LineStyle','--');
plot(t_vec,y_ex(t_vec),'LineWidth',1.5);
legend('FWD eu','BWD eu','Exact solution');
return