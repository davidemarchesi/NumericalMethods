function[xn,k,err_acc]=gs_sys(A,b,x0,tol,nmax)
% Gauss-Seidel resolutive method for linear systems:
% 
% Inputs:
% A = System matrix;
% b = known terms;
% x0 = starting vector (given by the initial conditions);
% tol = tolerance on the normalized residual;
% nmax = maximum number of iterations;
%
% Outputs:
% xn = solution obtained;
% k = number of iterations;
% err_acc = a vector that has inside the error of the i-th iteration (which
% basically is the position of the error in the vector);

% variables initialization:
n = length(b);
xn = zeros(n,1);
k = 0;

% to use the G-S method a squared matrix 'nxn' invertible (so with 
% diagonal elements ~=0), and a vector x0 of length 'n' is needed,
% otherwise the method is not applicable

if ((size(A,1)~=n) | (size(A,2)~=n) | (length(x0)~=n))
  error('Incompatible dimensions of the matrix (not squared) or of the vector!');
elseif (prod(diag(A)) == 0)
  error('The matrix A is not invertible!')
else
    xv = x0;
    % 'r' is the residual of the method
    r = b - A*x0;
    % 'err' is the error of the method normalized
    err = norm(r)/norm(b);
    % initialization of the vector 'err_acc'
    err_acc = [];
    while (err>tol && k<nmax)
        for i = 1:n 
            % The difference with the jacobi method lays here, where we do
            % not have in the first multiplication the k-th term 'xv', but
            % the k+1-th term 'xn' itself
            xn(i) = (b(i) - A(i,1:i-1)*xn(1:i-1) - A(i,i+1:n)*xv(i+1:n)) / A(i,i);
        end
        k = k+1;
        r = b - A*xn;
        err = norm(r)/norm(b);
        err_acc=[err_acc, err];
        xv = xn;
    end
end
end


