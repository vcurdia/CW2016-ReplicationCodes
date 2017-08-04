function testsemidefpos(Sigma)

% testsemidefpos
%
% Tests whether matrix Sigma is semi-definite positive
%
% ..............................................................................
%
% Created: September 12, 2011 by Vasco Curdia
% Copyright (2011) by Vasco Curdia


[Ur,Dr] = eig((Sigma+Sigma')/2);
Dr = diag(Dr);
tol = eps(max(Dr)) * length(Dr);
idxD = (abs(Dr) > tol);
Dr = Dr(idxD);
negeig = sum(Dr<0); % number of negative eigenvalues
if negeig~=0
  fprintf('Warning: Previous Hessian is not positive semi definite!\n');
  negeig
  Dr
end

