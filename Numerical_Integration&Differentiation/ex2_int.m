% In this code there is another function example to see how the various
% integration algorithms work (To do this, will be watched both the values
% of the integrals and of the errors with different numbers of
% subintervals).
% Function
x = linspace (0, pi/2 , 1000);
f =@(x) cos(x) .* exp(sin(x)) ;
% Numbers of subintervals
vect = 2.^(0:4);
% Extremes
a = 0;
b = pi/2;
% Declaration of the plot variables
int_mp_vec = [];
int_trap_vec = [];
int_simpson_vec = [];
vec_err_mp = [];
vec_err_trap = [];
vec_err_simpson = [];
% Exact integral value
Iexact = exp(1)-1;
% Calculation of the integral values and errors with a different number of
% subintervals k
for k = vect %(cycle which iterates at every k step)
    int_mp = mediump_int(a,b,k,f);
    int_trap = trap_int(a,b,k,f);
    int_simpson = simpson_int(a,b,k,f);
    
    int_mp_vec = [int_mp_vec; int_mp];
    err_mp = abs(int_mp - Iexact);
    vec_err_mp = [vec_err_mp; err_mp];
    
    int_trap_vec = [int_trap_vec; int_trap];
    err_trap = abs(int_trap - Iexact);
    vec_err_trap = [vec_err_trap; err_trap];
    
    int_simpson_vec = [int_simpson_vec; int_simpson];
    err_simpson = abs(int_simpson - Iexact);
    vec_err_simpson = [vec_err_simpson; err_simpson];
end

figure()
hold on;
xlabel('Number of subintervals');
ylabel('Integral values');
title('Num. Subint. vs Int. values');
plot(vect,int_mp_vec,vect,int_trap_vec,vect,int_simpson_vec,'LineWidth',1.4);
legend('MP int','Trap. int','Simpson int');
hold off;
figure()
xlabel('Number of subintervals');
ylabel('Absolute error');
title('Num. Subint. vs Abs. Error');
loglog(vect,vec_err_mp,vect,vec_err_trap,vect,vec_err_simpson,'LineWidth',1.4);
legend('MP err','Trap. err','Simpson err');
hold off;