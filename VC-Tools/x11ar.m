function ysa=x11ar(y,freq)

% Conversion of GAUSS code for seasonal adjustment
% GAUSS original from M. Watson
% X11AR.prc, mww, 3/29/00
% Carry out seasonal adjustment using X11 AR
% This applies the linear version of X11 to the
% series y, after it has been padded out with forecasts
% and backcasts constructed from an estimated AR model
% The current AR model is an AR(3) with lags 12 and 13
% added in the monthly case and an AR(4) in the quarterly case. 
% This can be changed by the vector ARVEC below.
%
% ..............................................................................
%
% Created: October29, 2002 by Vasco Curdia
% Updated: July 26, 2011 by Vasco Curdia
% 
% Copyright 2002-2011 by Vasco Curdia

if isnan(y)
    error('Y contains missing values in X11AR. proc not implemented for missing values');
end

% Pad series
if freq=='m'
    n=84;
    arvec=[1;2;3;12;13];
elseif freq=='q'
    n=28;
    arvec=[1;2;3;4];
else
    error('Data frequency must be either monthly (m) or quarterly (q)');
end

ypad=padar(y,n,arvec);

% X11 filter
x11=x11filt(freq);

% Construct Seasonally Adjusted Series
ysa=zeros(size(y,1),1);
for i=1:size(y,1)
    ysa(i)=ypad(i:i+2*n)'*x11;
end
