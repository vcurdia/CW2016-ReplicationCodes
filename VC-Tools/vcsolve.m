function varargout = vcsolve(sys,varargin)

% This small function is intended at simplifying the call for solve by
% simply calling on the system to be solved and, possibly, the unknowns
% instead of calling each of the equations of the system
%
% ..............................................................................
%
% Created: August 16, 2004 by Vasco Curdia
% Updated: July 26, 2011 by Vasco Curdia
% 
% Copyright 2004-2011 by Vasco Curdia

%--------------------------------------------------------------------------

comstr = 'solve(';
for j=1: max(size(sys))
    comstr = [comstr,'sys(',num2str(j),'),'];
end
for j=1:nargin-1
    eval(['syms ',char(varargin{j})])
    comstr = [comstr,char(varargin{j}),','];
end
comstr = [comstr(1:end-1),')'];
solvesol = eval(comstr);
varfields = fieldnames(solvesol);
nvar=size(varfields,1);
if nargout<=1
    varargout{1} = solvesol;
elseif nargout==nvar
	for j=1:nvar
        varargout{j} = getfield(solvesol,varfields{j});
	end
else
    error([int2str(nvar) ' variables does not match ' ...
             int2str(nargout) ' outputs.'])
end
