function Sigma = vccovdefpos(Sigma,tol,epsilon)

% vccovdefpos
%
% ..............................................................................
%
% Created: October 29, 2002 by Vasco Curdia
% Updated: July 26, 2011 by Vasco Curdia
% 
% Copyright 2002-2011 by Vasco Curdia

% if ~exist('tol','var'), tol = 10*eps(full(max(abs(diag(Sigma))))); end
if ~exist('tol','var'), tol = 1e-10; end
if ~exist('epsilon','var'), epsilon = 1e-16; end

% NumPrecision = 1e-12;
% Epsilon = 1e-10;
% [uu,dd,vv]=svd(ptT);
[uu,dd]=eig((Sigma+Sigma')/2);
% dd=round(dd/NumPrecision)*NumPrecision;
dd = diag(dd);
dd(abs(dd)<tol)=epsilon;
% dd(dd<tol)=tol*10;
dd = diag(dd);
Sigma = uu*dd*uu';
Sigma = (Sigma+Sigma')/2;
% [T,num] = cholcov(Sigma,0)

