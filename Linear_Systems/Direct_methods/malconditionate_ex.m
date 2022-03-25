% The following will be an example of the efficacy on the resolution of the
% LU factorization with a malconditionate matrix.
% We will use as an example the Hilbert matrix (a classic example 
% of malconditionate matrix), for which MATLAB has already implemented 
% a function [press help to see how it is built].
% Obv. the system considered will be Hx=b

% We start the example taking a matrix dimension of 5
n = 5;
% Here the command to create the Hilbert matrix
H=hilb(n);
% Initialization of the variables
L=zeros(n);
U=zeros(n);
P=zeros(n);
% Here the LU factorization of the matrix
[L,U,P]=lu(H);

% Firstly we calculate the vector b assuming we know the uknown vector x,
% which will be called xex (exactly-x), that will be xex=[1;...;1]
xex=ones(n,1);
b=H*xex;
% Initialization of the variables x and acc
x = [];
acc = 0;
% Now we will solve the system with the function 'linear_sys_LU'
[x,acc]=linear_sys_LU(H,b);
% The conditioning number will be given by acc, and the relative error
% will be:
rel_err=0;
% (The difference between xex and x)
diff=xex-x;
rel_err= norm(diff,2)/norm(xex,2);

fprintf('The conditioning number is %d \n', acc);
fprintf('The relative error is %d \n', rel_err);
%%
% Now we will repeat the process for n=2,...,15
% Initialization of the variables:
% vcond= vector with the conditioning number of every iteration
vcond=[];
% vdiff= vector with the norm of the difference of xex and x, of every
%        iteration
vdiff=[];
% vrel_err= vector with the relative error of every iteration
vrel_err=[];
% rel_err= relative error of the i-th iteration
rel_err=0;
% diff= difference of every iteration, between xex and x
diff=0;
% norm_diff= norm of diff
norm_diff=0;
% nit= vector with the iteration number
nit=[];
for i=2:16
    
    H=hilb(i);
    L=zeros(i);
    U=zeros(i);
    P=zeros(i);

    [L,U,P]=lu(H);

    xex=ones(i,1);
    b=H*xex;

    x = [];
    acc = 0;

    [x,acc]=linear_sys_LU(H,b);
    
    vcond=[vcond,acc];

    diff=xex-x;
    norm_diff=norm(diff,2);
    rel_err=norm_diff/norm(xex,2);
    vrel_err=[vrel_err,rel_err];
    
    nit=[nit,i];
end

% Now there's a plot in semilog scale of the conditioning number (blue),
% and of the relative error (red), in function of the n-th iteration
semilogy( nit,vcond ,'o', 'LineWidth', 2);
hold on;
semilogy(nit,vrel_err,'*','Linewidth',2);
xlabel('Iteration number');
legend('conditioning number','relative error');
