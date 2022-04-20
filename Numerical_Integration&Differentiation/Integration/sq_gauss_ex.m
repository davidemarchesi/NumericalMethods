% This code has the aim to analize the Gauss Squaring formulas.
% Let's consider a given function :
f =@(x) x.^6 +1;
% Analytically the exact integral in the interval [0,1] can be easily
% calculated, and the result is :
I_ex = 1/7 + 1;
% Now we want to analize how the error obtained with the calculation of
% this exact integral changes with different numbers of nodes using the
% implemented formula 'gauss_sq_I.m' up to a number of m=5 nodes.
%
% Interval definition :
a = 0;
b = 1;
% Declaration of the variables :
m = [1; 2; 3; 4; 5];
I_err = 0;
I_err_vec = [];

for k=m
    I = gauss_sq_I(a, b, f, m);
    I_err = abs(I_ex-I);
    I_err_vec = [I_err_vec ; I_err];
end

figure()
xlabel('Polynomial grade');
ylabel('Absolute error');
hold on;
title('Error path in relation with the number of nodes');
plot(m,I_err_vec,'LineWidth',1.4);
legend('m=1','m=2','m=3','m=4','m=5','Location','eastoutside');
return