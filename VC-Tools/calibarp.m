function [tMax,HalfLife]=calibarp(varargin)

% calibarp
%
% ..............................................................................
%
% Created: October 29, 2002 by Vasco Curdia
% Updated: July 26, 2011 by Vasco Curdia
% 
% Copyright 2002-2011 by Vasco Curdia

%% options
nSteps = 20;
ShockSize = 1;

%% evaluate inputs
p = nargin;
Phi = [varargin{:}; eye(p-1) zeros(p-1,1)];
Omega = [ShockSize;zeros(p-1,1)];
idx = eye(p,1);

%% simulate
z = zeros(p,nSteps);
z(:,1) = Omega; 
for t=2:nSteps
    z(:,t) = Phi*z(:,t-1);
end
y = idx'*z;

%% Plot
plot(1:nSteps,y)
xlim([1,nSteps])

%% Get max and halflife
[yMax,tMax] = max(y);
HalfLife = tMax+find(y(tMax+1:end)<=yMax/2,1,'first');

%% show results:
fprintf('\nAR(%.0f):',p)
fprintf('\n---------')
fprintf('\n     yMax: %.5f',yMax)
fprintf('\n     tMax: %.0f',tMax)
fprintf('\nHalf-life: %.0f\n\n',HalfLife)


