function [X,U] = Dir_Neuelliptic(mu,beta,c,fun,a,b,gN,h)

% [X,U] = Dir_Neuelliptic(mu,beta,c,fun,a,b,gN,h)
%
% This function solves the problem at limits (Omogeneus Dirichlet problem 
% on the left extreme, and Neumann non-omogeneus problem on the right extreme)
% with constant coefficients mu, beta, c. Here follows an example of a
% function that works with this script :
% 
%  -mu u'' + beta u' + cu = f(x), in (a,b),
%  u(a) = 0,
%  mu*u'(b) = g_N.
% 
% For the resolution is used the FEM (linear).
%
% INPUT:
% mu   = diffusion coefficient;
% beta = transport coefficient;
% c    = reaction coefficient;
% fun  = known term function handle;
% a    = left extreme on the 1D domain;
% b    = right extreme on the 1D domain;
% gN   = Neumann data on the right extreme;
% h    = discretization step.
% 
% OUTPUT:
% 
% X = vector which contains the grid nodes;
% U = vector which contains the approximated solution of the grid

% Border nodes
x0 = a;
xL = b;

% Internal nodes
xin = [x0+h :h :xL-h]';
N = length(xin);

% Stiffness Matrix (A)
d0  = ones(N,1);
d1  = ones(N-1,1);

Aint = mu/h  .* (2*diag(d0) - diag(d1,1) - diag(d1,-1))+ ...
    beta  .* (1/2*diag(d1,1) - 1/2*diag(d1,-1))+ ...
    c.*h   .* (2/3*diag(d0) + 1/6*diag(d1,1) + 1/6*diag(d1,-1));
% Definition of the matrix for the Neumann problem
A = zeros(N+1,N+1);
A(1:N,1:N) = Aint;
A(N+1,N) = mu/h * (-1) + beta * (-1/2) + c * h * (1/6);
A(N,N+1) = mu/h * (-1) + beta * ( 1/2) + c * h * (1/6);
A(N+1,N+1) = mu/h + beta * (1/2) + c*h * (1/3);

% Known term (f)
fint = h*fun(xin);
f = zeros(N+1,1);
f(1:N,1) = fint;
f(N+1,1) = h/2*fun(xL)+gN;

% Check of the symmetry of the matrix A (granted if beta=0 : NO TRANSPORT)
if A == A'
   fprintf('A is symmetric \n')
else
   fprintf('A is not symmetric \n')
end

% Check of the positivity of the matrix A 
% (granted if mu>0 & c>=0 : COERCIVITY OF a(u,v))
if min(eig(A))>0
    fprintf('A is positive definite \n')
else 
    fprintf('A is not positive definite \n')
end

% Resolution of the linear system (with the \ MATLAB function)
u = A\f;

% Addition of the borders for the grid nodes and for the resolutions
X = [x0;xin;xL];
U = [0;u];

end