function[x]=bksub(U,y)
% Backward substitution method for the resolution of the Upper triangular
% matrix from a square matrix A of a linear system. This function is also
% implemented to be used with the function 'linear_sys_LU'
n = length(U);
x = zeros(n,1);

x(n)=y(n)/U(n,n);

for  i = n-1:-1:1
    x(i) = (y(i) - ( U(i,i+1:n)*y(i+1:n)))/U(i,i);
end
end
