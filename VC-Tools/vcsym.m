function vcsym(varargin)

% vcsym
%
% Creates symbolic variable(s) in caller workspace.
%
% Usage:
%   vcsym(<VarName1>,...,<VarNameN>)
%
% Inputs:
%   <VarName> (character array)
%   Name of the symbolic variable to create.
%
% ...........................................................................
% 
% Created: September 29, 2016 by Vasco Curdia
% 
% Copyright 2016-2017 by Vasco Curdia

%% -------------------------------------------------------------------

for j=1:nargin
    assignin('caller',varargin{j},sym(varargin{j}))
end

%% -------------------------------------------------------------------

