function output = bwd (f,x,h)
% A simple script for the 'Backward' method [BW] (Finite difference method)
output = (f(x)-f(x-h))./h ;
end