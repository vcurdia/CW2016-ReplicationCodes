function vcdispvar(v,fid)

% vcdispvar
%
% Displays the names and values of a list of variables
%
% Usage:
%   vcDispVar(v)
%   vcDispVar(v,fid)
%
% with:
%   v: list of variable names to show. If an individual entry is empty,
%      using either '' or [], then an empty line is shown at that position.
%      The code will evaluate the variable name on the caller workspace to
%      extract the value to display.
%   fid: file identifier if displaying results into a text file.
%
% .........................................................................
%
% Created: December 31, 2012 by Vasco Curdia
%
% Copyright 2012 by Vasco Curdia

if ~exist('fid','var'),fid = 1;end
maxlength = max(cellfun('length',v));
for j=1:length(v);
  if isempty(v{j})
    fprintf(fid,'\n');
  else
    fprintf(fid,['%',int2str(maxlength),'s: %7.3f\n'],v{j},evalin('caller',v{j}));
  end
end
