function x=roundd(x,d)
% rounds the number x to the d-th decimal case
%
% ..............................................................................
%
% Created: June 17, 2004 by Vasco Curdia
% Updated: July 26, 2011 by Vasco Curdia
% 
% Copyright 2004-2011 by Vasco Curdia

x = round(x*10^d)*10^(-d);