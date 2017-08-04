% parallelgain
% 
% Computes the gains from using parallel computing, based on the number of
% parallel machines
%
% The idea is that by computing all possibilities of n successive
% iterations implies 2^n-1 computations, but which are shared by p
% processors, hence the number of operations per machine is only (2^n-1)/p.
%
% Computing iteratively implies n computations by a single machine.
%
% Therefore for any given p, the gain from parallel computing is given by 
% n-(2^n-1)/p
% 
% The optimal number of successive iterations to feed parallel computing at
% a time is given by maximizin the gain w.r.t. n, conditional on p, which
% yields the FOC: nOpt(p)=(log(p)-log(log(2)))/log(2).
%
% nMax(p) is the highest number of successive draws that still yields a
% positive gain.
%
% ..............................................................................
% 
% Created: August 2, 2007 by Vasco Curdia
% Updated: July 26, 2011 by Vasco Curdia
% 
% Copyright 2007-2011 by Vasco Curdia

% -------------------------------------------------------------------------

clear all
Np = 32;
for p=1:Np
    n=1:p;
    gain=n-(2.^n-1)/p;
    nMax(p)=max(n(gain>=0));
    nMaxGgain(p)=gain(nMax(p));
    nOpt(p)=(log(p)-log(log(2)))/log(2);
    nOptR(p)=max(round(nOpt(p)),1);
    nOptRGain(p)=gain(nOptR(p));
end
[1:Np;nOpt;nOptR;nOptRGain;nMax;nMaxGgain]'
