function t = addtime(t,tDelta)

% addtime
%
% adds tdelta periods to the date t.
%
% Usage:
%   addtime(t)
%   addtime(t,tdelta)
%
% Inputs:
%
%   t (string)
%   Time period. 
%
%   tDelta (string)
%   Number of time periods to add. Default: 1
%
% Output:
%
%   t (string)
%   New time.
%
% Convention used:
% Dates in format of '####q#' or '####m##' for quarterly and monthly data,
% respectively. For monthly frequency, for single digit periods they can be show
% up in TimeStart and TimeEnd as either # or ## (e.g. 1 or 01) but the index
% will be in the format ##.
%
% See also
% createTimeIdx
%
% ..............................................................................
% 
% Created: March 6, 2014 by Vasco Curdia
% 
% Copyright 2014-2016 by Vasco Curdia

%% -----------------------------------------------------------------------------

%% default
if nargin<2
    tDelta = 1;
end

%% Find frequency
if ismember('q',t)
    nPeriods = 4;
    PerStr = 'q';
    nPerStr = 1;
elseif ismember('m',t)
    nPeriods = 12;
    PerStr = 'm';
    nPerStr = 2;
else
    error('Frequency could not be detected.')
end

%% Generate new time
StartYear = eval(t(1:4));
StartPer = eval(t(6:end));
% NewPer = StartPer + tDelta;
% DeltaYear = floor(NewPer/nPeriods);
% NewYear = StartYear + DeltaYear;
% NewPer = NewPer-DeltaYear*nPeriods;
tStart = StartYear+(StartPer-1)/nPeriods;
tNew = tStart+tDelta/nPeriods;
NewYear = floor(tNew);
NewPer = 1+(tNew-NewYear)*nPeriods;
t = sprintf(['%04.0f%s%0',int2str(nPerStr),'.0f'],NewYear,PerStr,NewPer);

%% -----------------------------------------------------------------------------
