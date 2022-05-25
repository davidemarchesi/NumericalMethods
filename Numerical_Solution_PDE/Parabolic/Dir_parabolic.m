function [X,U,T,dt_max] = Dir_parabolic(mu,beta,c,fun,g,a,b,t0,t_max,h,tau,theta)

% [X,U,T,dt_max] = Dir_parabolic(mu,beta,c,fun,g,a,b,t0,t_max,h,tau,theta)
%
% Solves the omogeneus Cauchy-Dirichlet problem for the diffusion-transport-reaction
% equation with constant coefficients mu,beta,c. Here an example :
%
%  u_t -mu u_xx + beta u_x + cu = f(x,t), for x in [a,b], t in (0,t_max]
%  u(a,t) = u(b,t) = 0,                                   t in (0,t_max]
%  u(x,0) = g(x),                         for x in [a,b].
%
% Is used the FEM (linear) in space, and the theta-method in time.
% 
% INPUT:
% mu    = diffusion coefficient;
% beta  = transport coefficient;
% c     = reaction coefficient;
% fun   = known term function handle (depending on x and t);
% g     = initial condition function handle (depending on x);
% a     = left extreme on the 1D domain;
% b     = right domain on the 1D domain;
% t0    = initial time;
% t_max = final time;
% h     = spatial discretization step;
% tau   = temporal discretization step;
% theta = parameter to select the numerical method temporal advancement.
%         It can assume values in [0,1];
%
% OUTPUT:
% X = Vector which contains the grid nodes;
% U = matrix N_nodes x N_instants which contains the approximated solution
%     on the grid nodes for every temporal instant (for every instant the solution
%     on the grid nodes is a column vector);
% T = vectors with the temporal instants used for the solutions;
% dt_max = maximum temporal discretization step to grant that the condition of 
%          absolute stability is satisfied when using the Eulero forward method.

% Border nodes
x0 = a;
xL = b;

% Internal nodes
xin = [x0+h:h:xL-h]';
N = length(xin);

% Mass matrix (M)
d0  = ones(N,1);
d1  = ones(N-1,1);

M = h * (2/3*diag(d0) + 1/6*diag(d1,1) + 1/6*diag(d1,-1));

% Stiffness matrix (A)
A = mu/h  * (2*diag(d0) - diag(d1,1) - diag(d1,-1))+ ...
    beta  * (1/2*diag(d1,1) - 1/2*diag(d1,-1))+ ...
    c*h   * (2/3*diag(d0) + 1/6*diag(d1,1) + 1/6*diag(d1,-1));

% lambdamax=max(abs(eig(-M\A))) 

% Check of the A symmetry (granted if beta=0 , NO TRANSPORT)
if A == A'
   fprintf('A is symmetryc \n')
else
   fprintf('A is not symmetryc \n')
end

% Check of the A positivity (granted if mu > 0 & c >= 0 , COERCIVITY OF a(u,v))
if min(eig(A))>0
    fprintf('A is positive definite \n')
else 
    fprintf('A is not positive definite \n')
end

% dt_max to verify the condition of absolute stability (in case of explicit method)
dt_max = 2/max(abs(eig(A,-M)));
 
% matrix C1 = (1/dt)*M + theta*A
C1 = (1/tau)*M + theta*A;

% matrix C2 = (1/dt)*M + (theta-1)*A = C1 - A
C2 = C1 - A;

% Initial condition u_0 = G
G = g(xin);

% Vector with the temporal instants with which is solved the problem
T = [t0:tau:t_max];

% Initialization of the matrix that which will contain the approximated 
% solution of the internal notes xin
N_instants = length(T);
Uin = zeros(N,N_instants);

% Initial conditions imposition
Uin(:,1) = G;

% Iterative cycle that calculates : C1 * u_n = rhs, with 
% rhs = C2 * u_(n-1) + theta*F_n + (1-theta)*F_(n-1) in the internal nodes xin

% Resolution of the linear systems (in time) using the function "\" of MATLAB

for it = 2:N_instants
    
    u_old=Uin(:,it-1);
    
    t_old = T(it-1);
    t = T(it);
    
    F_old = h*fun(xin,t_old);
    F = h*fun(xin,t);
    
    rhs = C2*u_old + theta*F + (1-theta)*F_old;
    
    Uin(:,it) = C1\rhs;
    
end

% Addition of the borders for the grid nodes and for the solutions
X = [x0;xin;xL];
U = [zeros(1,N_instants);Uin;zeros(1,N_instants)];

end