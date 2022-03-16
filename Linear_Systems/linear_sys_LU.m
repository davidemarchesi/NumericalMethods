function[x,acc]=linear_sys_LU(A,b)
% Linear systems resolution with the LU factorization.
% The A matrix is factorized in LU matrices, and if necessary also using
% the Pivoting matrix (Ax=b --> PAx=LUx=Pb), and then the system is solved.
% This function NEEDS BOTH THE FUNCTIONS 'bksub' and 'fwsub' to run.
% Input:
% A = square matrix (n,n);
% b = vector with dimension (n,1);
% 
% Output:
% L = the lower triangular matrix of A;
% U = the upper triangular matrix of A;
% P = the pivoting matrix (If necessary, otherwise is not printed 
%     [it would be equal to the identity matrix]);
% x = the solution vector (n,1);
% acc = the conditioning of the matrix A, which gives us important
%       informations about the accuracy of the solution we have;

% First is good to see if the matrix is squared, in this way 
% if there was an error in the inserction of it, we do not even start the
% function
[n,m]= size (A);
if n~=m
    error('The matrix A IS NOT SQUARED! \n');
else
    
    % A factorization function is already implemented in MATLAB
    [L,U,P]=lu(A);
    % The identity matrix with dimension (n,n)
    I=eye(n);

    if P==I
        % here the case with the pivoting matrix being equal to the
        % identical one
        y=fwsub(L,b);
        x=bksub(U,y);

        fprintf('The Pivoting matrix is equal to the identity matrix I = eye (%d) \n', n);
        fprintf('The result is x = \n')
        disp(x);

    else
        % If the pivoting matrix is not identical, we have to insert it
        % into the 'fwsub' function
        Pb=P*b;
        y=fwsub(L,Pb);
        x=bksub(U,y);

        fprintf('The pivoting matrix is not the Identity matrix \n');
        fprintf('The result is x = \n');
        disp(x);
    end
    acc = cond(A);
    fprintf('\n The solution is related with a conditioning number of the matrix A of : %d \n', acc);
end
end
   
