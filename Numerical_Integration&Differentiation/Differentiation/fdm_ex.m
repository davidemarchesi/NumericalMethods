% This is a simple script to evidence the differences and the
% characteristics of different algorithm for the 'Finite Differences
% method'. Here the function that will be considered for the analysis :
f = @(x) exp(-x);
df = @(x) -exp(-x);
x_0 = 0.25;
df_x_0 = df(x_0);
%
% Now we want to see the esteem of the local derivative in the point x_0
% with each of the three different methods that can be used : BW , FW and
% CM. To do this we use the three functions implemented:
% (the value of 'h' must be set)
h = 0.001;
fw_dx = fwd(f,x_0,h);
bw_dx = bwd(f,x_0,h);
cm_dx = cmd(f,x_0,h);
fprintf('\n The value of the local derivate using the 3 different methods is respectively :\n');
fprintf(' %f for the FW; \n %f for the BW; \n %f for the CM. \n', fw_dx, bw_dx, cm_dx);
%%
% Now we want to verify the convergence grade for each algorithm variating
% the discretization interval.
%
% Variable declaration
err_fwd = 0;
err_bwd = 0;
err_cmd = 0;
err_fwd_vec = [];
err_bwd_vec = [];
err_cmd_vec = [];
h_vec = [0.2; 0.1; 0.05; 0.025; 0.0125];
fw_dx = 0;
bw_dx = 0;
cm_dx = 0;
for h = h_vec
    fw_dx = fwd(f,x_0,h);
    err_fwd = abs(df_x_0 - fw_dx);
    err_fwd_vec = [err_fwd_vec; err_fwd];

    bw_dx = bwd(f,x_0,h);
    err_bwd = abs(df_x_0 - bw_dx);
    err_bwd_vec = [err_bwd_vec; err_bwd];

    cm_dx = cmd(f,x_0,h);
    err_cmd = abs(df_x_0 - cm_dx);
    err_cmd_vec = [err_cmd_vec; err_cmd];

end

% To make an esteem of the convergence is sufficient to watch the pendences
% of the error lines in a logarithmic scale : the more is the slope, the
% quicker they converge. This will be done for each different algorithm 
% variating the step length 'h' :
figure()
loglog(h,err_fwd_vec,'LineWidth',1.4,'Color','r');
hold on
loglog(h,err_bwd_vec,'LineWidth',1.4,'Color','b');
% Line to compare with the convergence order
loglog(h,h,'LineWidth',1.5,'Color','g',LineStyle='-');
loglog(h,err_cmd_vec,'LineWidth',1.4,'Color','k');
% Squared function to compare with the convergence order
loglog(h,h.^2,'LineWidth',1.5,'Color','k',LineStyle='-.');
grid on;
xlabel('Step length');
ylabel('Absolute error');
title('Convergence and Error of different methods');
legend('FW','BW','1gr','CM','2gr','Location','eastoutside');
% Where is clear that the first two methods, the BW and the FW have a first
% order convergence, while the CM has a second order convergence.
%%
% Implementing the same script we have above for smaller steps 'h' some
% interesting observation can be made
err_fwdn = 0;
err_bwdn = 0;
err_cmdn = 0;
err_fwd_vecn = [];
err_bwd_vecn = [];
err_cmd_vecn = [];
vec1 = 0:50;
vec2 = 2.^(vec1);
h_vecn = 0.2./vec2;
fw_dxn = 0;
bw_dxn = 0;
cm_dxn = 0;
for h = h_vecn
    fw_dxn = fwd(f,x_0,h);
    err_fwdn = abs(df_x_0 - fw_dxn);
    err_fwd_vecn = [err_fwd_vecn; err_fwdn];

    bw_dxn = bwd(f,x_0,h);
    err_bwdn = abs(df_x_0 - bw_dxn);
    err_bwd_vecn = [err_bwd_vecn; err_bwdn];

    cm_dxn = cmd(f,x_0,h);
    err_cmdn = abs(df_x_0 - cm_dxn);
    err_cmd_vecn = [err_cmd_vecn; err_cmdn];

end

figure()
loglog(h_vecn,err_fwd_vecn,'LineWidth',1.4,'Color','r');
hold on
loglog(h_vecn,err_bwd_vecn,'LineWidth',1.4,'Color','b');
loglog(h_vecn,err_cmd_vecn,'LineWidth',1.4,'Color','k');
grid on;
xlabel('Step length');
ylabel('Absolute error');
title('Error of different methods');
legend('FW','BW','CM','Location','eastoutside');
% THIS IS VERY IMPORTANT
% As we can see from this last graphic, if for the Finite Difference method
% are used too small steps, then the prevalent error is no more a numerical
% one, it becomes the floating point annotation error and it increases
% really quickly; due to this, this type of method is used only for quick
% evaluation, for more precise results industrial codes uses different type
% of methods.
%%
% The same analysis made above can be also made for greater derivative
% orders, one could be for example the code 'cmd2.m'.