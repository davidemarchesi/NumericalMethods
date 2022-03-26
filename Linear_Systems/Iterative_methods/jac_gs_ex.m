% This will be a short exercise to use both the Jacobi and the Gauss-Seidel
% resolutive methods, and also to compare them.
%
% Lets form the following system Ax=b :
b=[7;4;5;5;5;4;7];
n=length(b);
C1=9;
C2=-3;
C3=1;
A= diag(ones(n,1))*C1 + diag(ones(n-1,1),1)*C2 + diag(ones(n-1,1),-1)*C2 ...
    + diag(ones(n-2,1),2)*C3 + diag(ones(n-2,1),-2)*C3;

% Now lets examinate the iterative matrices related to the Jacobi and G-S
% method :
%
% Here the Jacobi iterative matrix :
I = eye(n);
Dinv = diag(1./diag(A));
Bj = I-Dinv*A;
%
% Here the G-S iterative matrix :
D = diag(diag(A));
E = -(tril(A,-1));
F = -(triu(A,1));
P = D-E;
Pinv = P\I;
Bgs = Pinv*F;
%
% Now examinate their spectral radius :
eigenvalues_j = eig(Bj);
spec_r_j = max(abs(eigenvalues_j)); % jacobi spectral radius
eigenvalues_gs = eig(Bgs);
spec_r_gs = max(abs(eigenvalues_gs)); % gauss-seidel spectral radius
% The system can be solved with the method when the spectral radius 
% calculated is < 1

% Now lets see the solution of the system with both the methods and how the
% error and the iterations change from one to the other :
% (we assume the following inputs : )
tol = 10^-6;
x0 = zeros(n,1);
nmax = 1000;
% solution with Jacobi :
[xn_j,k_j,err_acc_j]=jacobi_sys(A,b,x0,tol,nmax);
% solution with G-S :
[xn_gs,k_gs,err_acc_gs]=gs_sys(A,b,x0,tol,nmax);
%
fprintf('The number of iterations of the G-S method is : %d \n',k_gs);
fprintf('\nThe number of iterations of the Jacobi method is : %d \n', k_j);

k_j_vec=[];
k_gs_vec=[];

for ii=1:k_j
    k_j_vec=[k_j_vec,ii];
end

for jj=1:k_gs
    k_gs_vec=[k_gs_vec,jj];
end

plot(k_j_vec,err_acc_j,'Marker','*','LineWidth',0.5,'Color','b','DisplayName','Jacobi');
hold on
plot(k_gs_vec,err_acc_gs,'Marker','o','LineWidth',0.5,'Color','r','DisplayName','G-S');
xlabel('Iteration num.');
ylabel('error at the iteration');
legend
