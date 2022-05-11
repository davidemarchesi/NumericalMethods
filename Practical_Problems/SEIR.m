% This following one will be the visualization of a practical ODE problem
% resolution. Will be approached the SEIR epidemiological model.
% For further informations about the problem refer to the pdf with a
% description of it.

% Assigned exstimated constant of the problem :
N0 = 10000;
t_0 = 0;
t_max = 100;
Tinf = 2;
Tinc = 5;

% Declaration of the initial conditions :
I0 = 0;
R0 = 0;
E0 = 1/N0;
S0 = 1-1/N0;

% Declaration of the three constant parameters which control the transition
% between one compartment and another :
beta = 1.75;
alpha = 0.2;
gamma = 0.5;

% Step selected for the solution of the ODE problem :
h = 0.01;

% Declaration of the starting values with a vector, and of the ODE system
% in the form of a function handle
y_0 = [S0; E0; I0; R0];
f =@(t,x) [-beta*x(1)*x(3); ...
            beta*x(1)*x(3)-alpha*x(2); ...
            alpha*x(2)-gamma*x(3); ...
            gamma*x(3)];
% BE AWARE that the i-th starting condition correspond to its i-th
% equation of the problem !

% Output declaration :
t_h = [];
u_h = [];

% Solution of the problem using the vectorial Eulero forward :
[t_h, u_h] = eulero_fwd_vec (f,t_0,t_max,y_0,h);

% Plot of the solution :
figure()
hold on
xlabel('Time (days)');
ylabel('Number of people');
plot(t_h,u_h(1,:),'LineWidth',1.4,'LineStyle','-','Color','g');
plot(t_h,u_h(2,:),'LineWidth',1.4,'LineStyle','--','Color','c');
plot(t_h,u_h(3,:),'LineWidth',1.4,'LineStyle','-.','Color','r');
plot(t_h,u_h(4,:),'LineWidth',1.4,'LineStyle',':','Color','b');
legend('Susceptible','Exposed','Infectious','Removed','Location','eastoutside');
%%
% The model could be analyzed with different values simply with a 'for'
% cycle, changing at every step a variable to see the behaviour of the
% model.
% (for example, but 'var' could be whatever parameter)
% var_vec = [...];
% for var = var_vec
%       f =@(t,x) [-var*beta*x(1)*x(3); ...
%                   var*beta*x(1)*x(3)-alpha*x(2); ...
%                   alpha*x(2)-gamma*x(3); ...
%                   gamma*x(3)];
%       .........
% end
%%
% An example could be the change of the problem with the introduction of a
% new variable 'rho' which could be, for example, a parameter which
% indicates the social distancing.
% The parameter will be introduced in the first two equations :
% dS(t) = -rho*beta*x(1)*x(3);
% dE(t) = rho*beta*x(1)*x(3)-alpha*x(2);
%

% (We assume to use the parameters announced above)
t_h_n = [];
u_h_n = [];

y_0 = [S0; E0; I0; R0];

rho_vec = [1; 0.75; 0.5; 0.25; 0];

for j = 1 : 5
    figure(j+1);
    hold on
    title('SEIR with rho constant');
    xlabel('Time (days)');
    ylabel('Number of people');
    f_v =@(t,x) [-rho_vec(j).*beta.*x(1).*x(3); ...
                  rho_vec(j).*beta.*x(1).*x(3)-alpha.*x(2); ...
                  alpha.*x(2)-gamma.*x(3); ...
                  gamma.*x(3)];
    [t_h_n, u_h_n] = eulero_fwd_vec (f_v,t_0,t_max,y_0,h);
    plot(t_h_n,u_h_n(1,:),'LineWidth',1.4,'LineStyle','-','Color','g');
    plot(t_h_n,u_h_n(2,:),'LineWidth',1.4,'LineStyle','--','Color','c');
    plot(t_h_n,u_h_n(3,:),'LineWidth',1.4,'LineStyle','-.','Color','r');
    plot(t_h_n,u_h_n(4,:),'LineWidth',1.4,'LineStyle',':','Color','b');
    legend('Susceptible','Exposed','Infectious','Removed','Location','eastoutside');
    hold off
end
return