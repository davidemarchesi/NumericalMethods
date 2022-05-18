function [X,U] = DirEllipt (mu,beta,c,fun,a,b,h)
%
% [X,U] = DirEllipt (mu,beta,c,fun,a,b,h)
% 
% Function which solves the Boundaries Problem (omogeneus Dirichlet, 1D) with
% constant coefficients 'mu,beta,c' :
%
% -mu u'' + beta u' + cu = f(x)   in [a,b];
% u(a) = u(b) = 0 ;
% 
% with the Finite Elements Method (Linear).
% 
% INPUTS :
% mu   = diffusion coefficient ;
% beta = transport coefficient ;
% c    = reaction coefficient ;
% fun  = function handle for the known term ;
% a    = left extreme in the 1D domain ;
% b    = right extreme in the 1D domain ;
% h    = Spatial discretization step .
% 
% OUTPUTS :
% X = vector which contains the grid points ;
% U = vector which contains the approximated solutions of the grid points ;

% As said above this function will consider only the particular case of
% omogeneus Dirichlet with a 1D dimension, and with constant coefficients.
% In the end the resolution will pass through the resolution of the linear
% system : Au=F . Where :
% - The matrix A will be formed by :
%   
%   A = (mu/h * K) + (beta * B) + (ch * C) ;
%               
%              |2  -1      |     |0 1/2        |     |2/3 1/6      |
%              |-1  ...    |     |1/2  ...     |     |1/6  ...     |
%   Where : K =|     ... -1|; B =|      ... 1/2|; C =|      ... 1/6| .
%              |      -1  2|     |        1/2 0|     |      1/6 2/3|
%                
% - The known term F , which in this case, using the trapezoid method will
%   be given by the values :
%
%   F_i = h * f(x_i) ;
%

% Boundary nodes :
x_0 = a ;
x_end = b ;

% Internal nodes definition :
x = x_0+h : h : x_end-h ;
x = x' ;
N = length(x) ;

% Definition of the matrix A :
K = mu/h * (diag(2*ones(N,1)) + diag(-1*ones(N-1,1),1) + diag(-1*ones(N-1,1),-1)) ;
B = beta/2 * (diag(ones(N-1,1),1) + diag(-1*ones(N-1,1),-1)) ;
C = c * h/6 * (diag(4*ones(N,1)) + diag(ones(N-1,1),1) + diag(ones(N-1,1),-1)) ;
A = K+B+C ;

% Definition of the known term :
f = h.*fun(x) ;

% Check of some conditions of the matrix A :
% Check over the symmetry of A (granted if beta = 0 --> NO TRANSPORT)
if A == A'
   fprintf('A is symmetric \n')
else
   fprintf('A is NOT symmetric \n')
end

% Check on the positivity of A (granted if mu > 0, c >= 0 (conditions to satisfy the coercivity of a(u,v)))
if min(eig(A))>0
    fprintf('A is positive definite \n')
else 
    fprintf('A is NOT positive definite \n')
end

% Resolution of the linear system :
u = A\f ;

% Definition of the outputs adding the Dirichlet condition :
U = [0; u; 0] ;
X = [x_0; x; x_end] ;

return