function TimeIdx = mergetimeidx(t1,t2)

% mergetimeidx
%
% Merges two time indeces created by TimeIdxCreate
%
% Usage:
%   mergetimeidx(t1,t2)
%
% Inputs:
%
%   t1 (cell)
%   First time index. 
%
%   t2 (cell)
%   Second time index.
%
% Output:
%
%   TimeIdx (cell array of strings)
%   Time index.
%
% Convention used: Dates in format of '####q#' or '####m##' for quarterly and 
% monthly data, respectively. For monthly frequency, for single digit periods 
% they can be show up in TimeStart and TimeEnd as either # or ## (e.g. 1 or 01)
% but the index will be in the format ##.
%
% See also
% createtimeidx
%
% .............................................................................
% 
% Created: June 25, 2014 by Vasco Curdia
% 
% Copyright 2014 by Vasco Curdia

%------------------------------------------------------------------------------

%% Match frequencies
if ~strcmp(t1{1}(5),t2{1}(5))
    error('Frequency of two indeces does not match.')
end

%% Find limiting dates
Bounds = [t1([1,end]), t2([1,end])];
Dates = zeros(1,4);
for j=1:4
    Dates(j) = eval(Bounds{j}(1:4)) + eval(Bounds{j}(6:end));
end
[~,idx] = sort(Dates);
Bounds = Bounds(idx([1,4]));

%% Create new index
TimeIdx = TimeIdxCreate(Bounds{:});

%------------------------------------------------------------------------------


