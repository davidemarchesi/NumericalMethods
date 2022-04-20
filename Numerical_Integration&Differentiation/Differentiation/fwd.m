function output = fwd (f,x,h)
% A simple script for the 'Forward' method [FW] (Finite difference method)
output = (f(x+h)-f(x))./h ;
end