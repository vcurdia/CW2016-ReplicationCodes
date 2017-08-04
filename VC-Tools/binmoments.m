function binmoments(m)
% computes the first m moments of a binomial distribution
% must have m>1
%
% ..............................................................................
%
% Created: April 16, 2004 by Vasco Curdia
% Updated: July 26, 2011 by Vasco Curdia
% 
% Copyright 2002-2011 by Vasco Curdia

syms n p t

G=(p*exp(t)+1-p)^n
moment=jacobian(G,t);
for j=2:m
    moment(j,1)=jacobian(moment(j-1),t);
end

moment = simplify(subs(moment,t,0))