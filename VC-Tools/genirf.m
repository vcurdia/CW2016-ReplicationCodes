function IRF = genirf(G0,G1,G2,T,Shocks)

% genirf
%
% Computes impulse responses to shocks. It assumes that the state space is given
% by:
%   z(t) = G0 + G1*z(t-1) + G2*e(t)
%
% Usage:
%    IRF = GenIRF(G0,G1,G2,T)
%    IRF = GenIRF(G0,G1,G2,T,Shocks)
%
% Inputs:
%
% G0,G1,G2
% matrices with the coefficients for the state space. In a true IRF the constant
% term is ignored, so G0 needs to be specified as a column of zeros. If G0 is
% specified as [], a column of zeros is used.
%
% T
% number of steps to compute
%
% Shocks (optional)
% Array with different combinations of e(t) shocks. If not specified, an identity
% matrix is assumed. If Shocks is specified as a 3 dimensional object, it is
% assumed that it has dimensions (ne x nShocks x TShocks) where ne is the number
% of exogenous variables, e(t), in the system, nShocks is the number of columns
% in Shocks (different combinations of e(t)), and TShocks is the number of
% periods over which there are shocks to the system. This allows the
% specification of sequences of shocks. Notice that if only one combination of
% shocks for multiple periods then need to set Shocks to be (ne x 1 x TShocks).
%
% Outputs:
%
% IRF 
% Array with dimensions (nVars x nSteps x nShocks) where nVars is the number of 
% variables in the system, nSteps, the number of steps computed for % the IRF 
% and nShocks the number of columns in Shocks. So IRF(1,2,3) corresponds to the 
% second period response of the first variable to the shocks in column 3.
%
% ..............................................................................
% 
% Created: December 17, 2010 by Vasco Curdia
% Updated: July 26, 2011 by Vasco Curdia
% 
% Copyright 2002-2011 by Vasco Curdia

%% -----------------------------------------------------------------------------

nVars = size(G1,1);
if isempty(G0), G0 = zeros(nVars,1); end
if ~exist('Shocks','var'),Shocks = eye(size(G2,2)); end
nShocks = size(Shocks,2);
TShocks = size(Shocks,3);
IRF = zeros(nVars,nShocks,T);
IRF(:,:,1) = G0+G2*Shocks(:,:,1);
for t=2:T
    if t<=TShocks
        IRF(:,:,t) = G0+G1*IRF(:,:,t-1)+G2*Shocks(:,:,t);
    else
        IRF(:,:,t) = G0+G1*IRF(:,:,t-1);
    end
end
IRF = permute(IRF,[1,3,2]);

%% -----------------------------------------------------------------------------
