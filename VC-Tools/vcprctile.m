function Y = vcprctile(X,p,dim)

% uses the prctile command but for any number of dimensions
%
% dim is the dimension along which we want to calculate the percentile
%
% ..............................................................................
%
% Created: May 23, 2005 by Vasco Curdia
% Updated: July 26, 2011 by Vasco Curdia
% 
% Copyright 2005-2011 by Vasco Curdia

%--------------------------------------------------------------------------

Xs = size(X);
Xn = ndims(X);
Xnidx = 1:Xn;
Xnidx1 = [dim,Xnidx(find(Xnidx~=dim))];
Xnidx2 = Xnidx(find(Xnidx~=dim));

Xp = permute(X,Xnidx1);
Xrs = reshape(Xp,[Xs(dim),prod(Xs(Xnidx2))]);

pXrs = prctile(Xrs,p);

Y = reshape(pXrs,[Xs(Xnidx2)]);
