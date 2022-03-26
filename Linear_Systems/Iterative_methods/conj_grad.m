function [xn,k,xdif,err] = conj_grad(A,b,xv,tol,nmax)
% Conjugate Gradient method for the resolution of linear systems
%
% Inputs :
% A = system matrix;
% b = known terms vector;
% xv = starting vector;
% nmax = maximum number of iterations;
% tol = tolerance;
%
% Outputs :
% xn = solution vector;
% k = number of iterations made;
% xdif = vector with the norm of every iteration advancement;
% err = vector with the relative errors on the residual;

% Firstly is important to control that we have in input a squared matrix A:
size1 = 0;
size2 = 0;
[size1,size2] = size(A);
if size1~=size2
    error('\nThe matrix "A" inserted is not squared ! \n');
end

% Initialization of the variables :
xn=[];
xdif=[];
err=[];
k=0;

bnrm2= norm(b);
if bnrm2==0.0
    bnrm2 = 1.0;
end

% definition of the residual and of the error :
r = b - A*xv;
error_norm = norm(r) / bnrm2;
err = [error_norm];

% If the error is already under the tolerance we don't even start the
% iterations and we stop here :
if error_norm < tol
    fprintf('The value "xv" inserted is already under the set tolerance! ');
    return
end

for k=1:nmax
    z = r;
    rho = r'*z;
    if k==1
        p = z;
    else
        beta = rho / rho_1;
        p = z + beta*p;
    end
    q = A*p;
    alpha = rho / (p'*q);
    xn = xv + alpha*p;
    r = r - alpha*q;
    error_norm = norm(r) / bnrm2;
    err = [err,error_norm];
    xdif = [xdif,norm(xn-xv)];
    xv = xn;
    if error_norm < tol 
       break
    end
    rho_1 = rho;
end

% See if the method has converged in the maximum number of iterations set :
if k==nmax
    fprintf('\nThe method does not converge in the maximum number of iterations %d', nmax);
else
    fprintf('\nThe method converges in %d iterations \n', k);
end
end