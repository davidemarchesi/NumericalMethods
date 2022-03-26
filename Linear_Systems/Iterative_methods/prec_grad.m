function [xn,k,err] = prec_grad(A,b,P,x0,tol,nmax)
%
% Preconditionate Gradient Method (Preconditionate Dynamic Richardson)
%
% Inputs :
% A = system matrix;
% b = known terms vector;
% P = preconditionator;
% x0 = starting vector;
% tol = tolerance;
% nmax = maximum number of iterations admitted;
%
% Outputs :
% xn = solution;
% k = number of iterations;
% err = normalized residual error;

% Firstly is important to control that we have in input a squared matrix A:
size1 = 0;
size2 = 0;
[size1,size2] = size(A);
if size1~=size2
    error('\nThe matrix "A" inserted is not squared ! \n');
end

% Initialization of the variables:
k = 0;
xn = x0;
err = [];
for k = 1:nmax
    r = b - A * xn;
    err_new = norm(r) / norm(b);
    err = [err, err_new];
    
    % Conditions verification :
    if err_new < tol
        fprintf('Richardson converges in %d iterations \n', k);
        return
    elseif k == nmax
        fprintf('Richardson does not converge in %d iterations \n', nmax)
        return
    end
    
    z = P\r;
%     if (nargin == 6)
        alpha = (z'*r) / (z'*A*z);  % alpha dynamic
%     end
    xn = xn + alpha * z;    
end
end