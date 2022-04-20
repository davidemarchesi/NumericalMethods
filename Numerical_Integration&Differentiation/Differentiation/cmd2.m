function output = cmd2 (f,x,h)
% A simple script for the 'Centered' method for a 2nd derivative [CM2] (Finite difference method)
output = (f(x+h)-2.*f(x)+f(x-h))./(h.^2) ;
end