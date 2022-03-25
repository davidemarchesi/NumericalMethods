% This will be a short exercise to use the jacobi method
% Firstly lets create a matrix of length 'n x n' s.t.
%%
% In order to visualize the type of matrix produced in the following code
% lines, run this paragraph:
n=5;
R1=1;
R2=-2;
A=diag(ones(n,1))*R2 + diag(ones(n-1,1),-1)*R1;
A(1,:)=1;
disp(A);
%%
n=100;
R1=1;
R2=-2;
A=diag(ones(n,1))*R2 + diag(ones(n-1,1),-1)*R1;
A(1,:)=1;

% The iterative Jacobi's matrix can be easily visualized as:
D= diag(diag(A));
I= zeros(n,n)+ diag(ones(n,1));
% obs: the Identical matrix could have been also produced with the command
% I=eye(n)
Dinv= diag(1./diag(A));
% and Bj:
Bj=I-(Dinv*A);

% The system Ax=b can be easily solved using the function 'jacobi_sys' with
% the following values (for example):

% outputs initialization:
k=0;
err_acc=[];
xn=[];
% inputs initialization:
b=ones(n,1);
b(1)=2;
x0=zeros(n,1);
tol=10^-6;
nmax=1000;
[xn,k,err_acc]=jacobi_sys(A,b,x0,tol,nmax);

% lets see how the error changes at every iteration:
kfin=k;
k_vec=[];
for ii=1:kfin
    k_vec = [k_vec, ii];
end
plot(k_vec,err_acc,'LineWidth',2);
xlabel('number of iterations');
ylabel('error at the correspondig iteration');


