function []=vcnormsim(mx,vx,bin)

% simulates the shape of a normal distribution given the mean and variance
%
% ..............................................................................
%
% Created: November 21, 2004 by Vasco Curdia
% Updated: July 26, 2011 by Vasco Curdia
% 
% Copyright 2004-2011 by Vasco Curdia

%--------------------------------------------------------------------------

x = normrnd(mx,sqrt(vx),1e+6,1);
figure
hist(x,bin)
fprintf('probability mass below zero: %f',normcdf(0,mx,sqrt(vx)))