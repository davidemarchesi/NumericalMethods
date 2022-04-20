% In this code section will be taken in consideration a fixed function f(x)
% to see the different behaviour of different Ingration algorithm
% 
% Definition of a fixed function
f=@(x) (x./(2*pi)).*sin(x);
% We will consider the interval [0,2pi]
x = linspace(0,2*pi,1000);
% Funtion plot
figure()
hold on;
plot(x,f(x),'LineWidth',1.5,'Color','b');
xlabel('x');
ylabel('y');
title('Function plot');
hold off;
%%
% Now we will use the three different algorithms, each one for a particular
% type of integration method
% trap_int(a,b,N,f) , for the Trapezium formula
% medium_int(a,b,N,f) , for the Medium point formula
% simpson_int(a,b,N,f) , for the Simpson formula
%
% Now will be considered a different number of subintervals for each
% iteration in order to see how the calculation changes with this number;
% In order to do this, knowing the exact integral of the function in the
% set interval, will be plotted the error values in a logarithmic scale.
%
% Number of the subinterval each iteration [0,2,4,8,16]
vect = 2.^(0:4);
% Definition of the extremes [a,b]
a = 0;
b = 2*pi;
% Definition of the plot variables
vec_err_mp = [];
vec_err_trap = [];
vec_err_simpson = [];
% Exact value of the function integral in the set interval
Iexact = -1;
for k = vect %(cycle which iterates at every k step)
    int_mp = mediump_int(a,b,k,f);
    int_trap = trap_int(a,b,k,f);
    int_simpson = simpson_int(a,b,k,f);
    
    err_mp = abs(int_mp - Iexact);
    vec_err_mp = [vec_err_mp; err_mp];
    
    err_trap = abs(int_trap - Iexact);
    vec_err_trap = [vec_err_trap; err_trap];

    err_simpson = abs(int_simpson - Iexact);
    vec_err_simpson = [vec_err_simpson; err_simpson];
end

figure();
title('Abs. Error vs Number of Subintervals');
xlabel('Number of Subintervals');
ylabel('Abs. Error value');
loglog(vect,vec_err_mp,'Color','b','LineWidth',1.4);
loglog(vect,vec_err_trap,'Color','g','LineWidth',1.4);
loglog(vect,vec_err_simpson,'Color','r','LineWidth',1.4);
legend('MP err','Trap. err','Simpson err');
hold off;