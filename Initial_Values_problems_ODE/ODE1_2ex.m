% In these script will be analysed the convergence property of both the
% Backward and the Forward Eulero's methods
%
% Let's consider the following Cauchy problem :
s = -5;
dy =@(y,t) s.*y ;
y_ex =@(t) exp(s*t) ;
t_0 = 0 ;
t_max = 10 ;
y_0 = 1;
t_vec = t_0:0.01:t_max;
% From the theory we know that the critic value of 'h' of the time steps
% length is :
h_critic = 2 /abs(s) ;
% Now we want to see what happens with the condition of criticity, with
% values greater, and with values lower
% 
% In order to do this we will simply sum or subtract small quantities
h_1 = h_critic + 0.05.*h_critic ;
h_2 = h_critic ;
h_3 = h_critic - 0.05.*h_critic ;
% Now let's observe what happens solving with these 3 values :
[t_h1,u_h1] = eulero_fwd(dy,t_0,t_max,y_0,h_1);
[t_h2,u_h2] = eulero_fwd(dy,t_0,t_max,y_0,h_2);
[t_h3,u_h3] = eulero_fwd(dy,t_0,t_max,y_0,h_3);
figure()
xlabel('time');
ylabel('functions');
title('Behaviour of the solution with a critical step "h", higher and lower');
hold on;
plot(t_h1,u_h1,'LineWidth',1.4,'LineStyle','--');
plot(t_h2,u_h2,'LineWidth',1.4,'LineStyle','-.');
plot(t_h3,u_h3,'LineWidth',1.4,'LineStyle',':');
plot(t_vec,y_ex(t_vec),'LineWidth',1.4,'LineStyle','-');
legend('over-critic','critic','sub-critic','exact solution');
hold off;
%%
% For Backward Eulero, the situation is different : it is a convergent
% method for all the 'h' taken; but using the fixed point ,ethod this is
% not true at all as to be convergent this other method needs that :
% |d/du (h*f(t,u))| < 1
% (Plot are not discussed here)