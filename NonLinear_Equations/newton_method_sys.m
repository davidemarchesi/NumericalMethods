function [ xvect , it ] = newton_method_sys ( x0 , nmax , tol , f , J )
% This function has the aim to use the Newton method to solve non-linear
% systems.
%
% Input parameters:
% x0 : starting point for the analysis;
% nmax : maximum number of iterations admitted;
% tol : the tolerance we want for the solution;
% f : which is the array containing all the non-linear equations of the
%     system we want to solve;
% J : which is the Jacobian of 'f';
%
% Output parameters:
% xvect : iterations results (the last is the solution);
% it : number of iterations;

% Ouput and variables inizialization
it=0;
xvect=[];
res=tol + 1;
xold=[];

while ( it<nmax && res>tol )
    xold = x0;
    if abs(det(J(x0)))< eps
        disp(['The cycle was interrupt due to the fact that the det. of the ' ...
            'Jacobian is zero \n']);
        it = it +1;
        break
    end
    % '\' is a function already implemented by MATLAB to solve linear
    % systems of equations
    x0 = x0 - J(x0)\f(x0);
    xvect = [xvect ; x0];
    % The residual is calculated using the norm of the difference between
    % the precedent value of the vector of the solutions and the precedent
    res = norm(xold - x0, "inf" );
    it = it +1;
end

if it==nmax
    disp(['The cycle was stopped because was  reached the maximum number of ' ...
        'iterations k = nmax = %d \n', it]);
else
    disp('The number of iterations made was : k = %d \n', it);
end

fprintf('The result is: \n');
xvect

return