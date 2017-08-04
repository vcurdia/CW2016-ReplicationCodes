function W=vcwiener(Ts,ns)

% generates a wienner process with T observations and ns times
% the output is a matrix Tsxns
%
% ..............................................................................
%
% Created: October 12, 2002 by Vasco Curdia
% Updated: July 26, 2011 by Vasco Curdia
% 
% Copyright 2002-2011 by Vasco Curdia

E=normrnd(0,1,Ts,ns);
W=cumsum(E)./sqrt(Ts);
