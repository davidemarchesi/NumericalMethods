function [succ, it] = fixpt(x_0, phi, nmax, tol)
% Fixed point method for the resolution of non-linear equations 
%
% Inputs :
% - x_0 : Starting Point;
% - phi : Fixed point function (inline or anonimus definite);
% - nmax : Maximum number of iterations;
% - tol : tolerance;
%
% Outputs :
% - succ : Vector containing all the iterations evaluated (the last is the
%                                                          solution);
% - it : iterations;
err   = 1 + tol;
it    = 0;
succ  = x_0;
xv    = x_0;
while (it < nmax && err > tol)
   xn    = phi(xv);
   err   = abs(xn - xv);
   succ = [succ; xn];
   it    = it + 1;
   xv    = xn;
end
%fprintf(' \n Numero di Iterazioni    : %d \n',it);
%fprintf(' Punto fisso calcolato   : %12.13f \n',succ(end));
end