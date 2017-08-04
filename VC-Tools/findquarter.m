function idx = findquarter(tid,q)

% findquarter
%
% Finds instances of quarter q in timeidx list tid. If not specified, first 
% quarter is assumed.
%
% See also:
% timeidx
% 
% ...........................................................................
% 
% Created: April 15, 2017 by Vasco Curdia
% 
% Copyright 2017 by Vasco Curdia

if nargin<2, q = 1; end

idx = find(cellfun(@(s)strcmp(s(end),int2str(q)),tid));
