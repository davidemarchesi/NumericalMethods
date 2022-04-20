function output = cmd (f,x,h)
% A simple script for the 'Centered' method [CM] (Finite difference method)
output = (f(x+h)-f(x-h))./(2*h) ;
end