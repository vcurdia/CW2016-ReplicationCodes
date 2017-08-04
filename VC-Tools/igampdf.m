function p=igampdf(x,a,b)
% evaluates the pdf of a inverse gamma distribution, using the scale factor (b)
% as matlab uses it for the gamma distribution
%
% ..............................................................................
%
% Created: March 10, 2004 by Vasco Curdia
% Updated: July 26, 2011 by Vasco Curdia
% 
% Copyright 2004-2011 by Vasco Curdia and Daria Finocchiaro

% -------------------------------------------------------------------------

p = zeros(size(x));
idx = x>0;
p(idx) = b^(-a)*x(idx).^(-(a+1)).*exp(-1./(b.*x(idx)))/gamma(a);