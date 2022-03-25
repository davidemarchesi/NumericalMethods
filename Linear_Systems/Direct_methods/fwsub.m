function[y]=fwsub(L,b)
% Forward substitution method for the resolution of the Lower triangular
% matrix from a square matrix A of a linear system. This function is also
% implemented to be used with the function 'linear_sys_LU'
n = length(L);
y = zeros(n,1);

y(1)=b(1)/L(1,1);

for i = 2:n
    y(i) = (b(i) - ( L(i,1:i-1)*y(1:i-1)))/L(i,i);
end
end
