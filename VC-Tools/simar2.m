function AR2 = simar2(rho1,rho2,nSteps,SmallNumber)

% simar2
%
% Simulate AR(2) IRF.
%
% Usage: 
%
%   AR2 = simar2(rho1,rho2)
%   ... = simar2(rho1,rho2,nSteps)
%   ... = simar2(rho1,rho2,nSteps,SmallNumber)
%
% If no arguments are given then a default example is used.
%
% ..............................................................................
%
% Created: January 5, 2014 by Vasco Curdia
%
% Copyright 2014 by Vasco Curdia

%% -----------------------------------------------------------------------------

%% Options

if nargin<1, rho1 = 1.5; end
if nargin<2, rho2 = -0.55; end
if ~exist('nSteps','var') || isempty(nSteps), nSteps = 100; end
if ~exist('SmallNumber','var') || isempty(SmallNumber), SmallNumber = 1e-2; end

%% Make IRF

tid = 1:nSteps;
x = zeros(1,nSteps);

x(1) = 1;
x(2) = rho1*x(1);
for t=3:nSteps
  x(t) = rho1*x(t-1) + rho2*x(t-2);
end

%% Plot
plot(tid,x)
axis tight

%% Store results
AR2.rho1 = rho1;
AR2.rho2 = rho2;
[AR2.Max,AR2.tMax] = max(x);
AR2.Halflife = find((x./AR2.Max)>=0.5,1,'last');
AR2.End = find((x./AR2.Max)>=SmallNumber,1,'last');

%% -----------------------------------------------------------------------------
