% This exercise is to visualize the gradient methods for linear systems
%
% We will assume the following data and the following system :
tol = 10^-5;
nmax = 10000;
n = 50;
x0 = zeros(n,1);
C1 = 0.2;
b = ones(n,1)*C1;
C2 = 4;
A = diag(ones(n,1))*C2 - diag(ones(n-1,1),1) - diag(ones(n-1,1),-1) ...
    - diag(ones(n-2,1),2) - diag(ones(n-2,1),-2);

% Now we will calculate the conditioning number without using the 'cond'
% command already implemented in MatLab
I = eye (n);
Ainv = A\I;
cond_num = norm(A) * norm(Ainv);
% (with the MatLab command it would have been : 'cond_num = cond(A)')
fprintf('\nThe conditioning number of the matrix "A" is : %d \n', cond_num);

% Below, the resolution of the linear system with the three different
% methods :
%
% Classical gradient method
[xn_dyn,k_dyn,xdif_dyn,err_dyn]=dyn_grad(A,b,x0,tol,nmax);
ii=0;
k_dyn_vec=[];
for ii=1:k_dyn+1
    k_dyn_vec=[k_dyn_vec,ii];
end
%
% Conjugate gradient method
[xn_conj,k_conj,xdif_conj,err_conj]=conj_grad(A,b,x0,tol,nmax);
jj=0;
k_conj_vec=[];
for jj=1:k_conj+1
    k_conj_vec=[k_conj_vec,jj];
end
%
% Preconditionate gradient method
% (Is defined the preconditionator P)
C3 = 2;
P = diag(ones(n,1))*C3 - diag(ones(n-1,1),1) - diag(ones(n-1,1),-1);
[xn_prec,k_prec,err_prec]=prec_grad(A,b,P,x0,tol,nmax);
tt=0;
k_prec_vec=[];
for tt=1:k_prec
    k_prec_vec=[k_prec_vec,tt];
end

fprintf('\nThe Dynamic gradient method, does %d iterations\n', k_dyn);
fprintf('\nThe Conjugate gradient method, does %d iterations\n', k_conj);
fprintf('\nThe Preconditionate gradient method, does %d iterations\n', k_prec);
% Here follows a graphic :
plot(k_dyn_vec,err_dyn,'LineWidth',1.2,'Color','b','DisplayName','Dynamic');
hold on
plot(k_conj_vec,err_conj,'LineWidth',1.2,'Marker','*','Color','r','DisplayName','Conjugate');
hold on
plot(k_prec_vec,err_prec,'LineWidth',1.2,'Marker','+','Color','g','DisplayName','Preconditionate');
xlabel('Iteration num.');
ylabel('Error at the iteration');
legend
