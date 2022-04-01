% This code has the aim to show the interpolation process in MATLAB using
% the functions already implemented in MATLAB  called 'polyfit' and
% 'polyval', using the Lagrange's Polynomial.
%%
% In order to do this we will consider the following function
f = @(x) x.*sin(x);
% We will consider the interval [-2,6]
x = linspace(-2,6);
% The plot of the function will be
figure();
plot(x,f(x),'LineWidth',1.75,'Color','k');
%%
% Now we are going to build the Interpolant Lagrange Polynomials with
% grades n=2,4,6 with equispaced nodes
hold on
title('Exact function vs Lagrange interpolating polynomial');
xlabel('x-values');
ylabel('y-values');
for n=[2,4,6]
    % To use the Lagrange polynomial of grade n, we have to take n+1 points
    xx = linspace(-2 , 6, n+1 );
    y = f(xx);
    % After we have created the points, the function polyfit gives us the
    % coefficients of a polynomial of grade 'n' passing in the points given
    % by 'xx' and 'y'
    K = polyfit(xx,y,n);
    % The function polyval calculates the values of the polynomial with the
    % coefficients 'K' in the points of the vector 'x'
    P = polyval(K,x);
    plot(x,P,'LineWidth',1.25);
    % This second plot indicates us the exact points extracted to form the
    % polynomials
    plot(xx,y,'o','LineWidth',1.85);
end
legend('Exact','P2','P2 points','P4','P4 points','P6','P6 points');
hold off
%%
% Here the plot of the errors path
figure();
title('Error in the interpolants path')
xlabel('x-values');
ylabel('Absolute error');
hold on
for n=[2,4,6]
    xx = linspace(-2 , 6, n+1 );
    y = f(xx);
    K = polyfit(xx,y,n);
    P = polyval(K,x);
    err = abs(f(x) - P);
    max_err = max(err);
    % In this plot we can see that the errors variates a lot, but that all
    % obviously are zero in the interpolating nodes
    plot(x,err,'LineWidth',1.25);
    fprintf('The maximum error of the Lagrange interpolating polynomial "P%d" is %f \n', n, max_err);
end
legend('P2','P4','P6');