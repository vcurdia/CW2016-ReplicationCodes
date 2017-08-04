function y=isvector(A)

% isvector
%
% ..............................................................................
%
% Created: May 23, 2005 by Vasco Curdia
% Updated: July 26, 2011 by Vasco Curdia
% 
% Copyright 2005-2011 by Vasco Curdia

%--------------------------------------------------------------------------

nd = ndims(A);

if nd==2
    [n1,n2]=size(A);
    y = (n1-1)*(n2-1);
    if y==0
        y = 1;
    else
        y = 0;
    end
else
    y = 0;
end
    