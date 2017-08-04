function [eq,vars,eq_jacob]=vceqjacob(eq,var_name)

% vceqjacob
%
% computes and presents the different coefficients in a linear symbolic
% expression
%
% sintax: [eq,vars,eq_jacob]=vceqjacob(eq,var_name)
%
% inputs:
%   - eq: the equation to be analyzed
%   - var_name: [optional] the symbolic variable name to be used as
%               normalizing factor for the coefficients
%
% outputs:
%   - eq: the equation, after normalization, if any
%   - vars: vector of symbolic variables present in the equation
%   - eq_jacob: the vector of coefficients
%
% ..............................................................................
%
% Created: August 30, 2006 by Vasco Curdia
% Updated: July 26, 2011 by Vasco Curdia
% 
% Copyright 2006-2011 by Vasco Curdia

%% ************************************************************************

vars = evalin('base',['[',findsym(eq),']']);
eq_jacob = double(jacobian(eq,vars));
if nargin == 2
    idx_var = double(jacobian(vars,var_name));
    if ~all(idx_var==0)
        idx_var = find(idx_var~=0);
        eq = eq/eq_jacob(idx_var);
        eq_jacob = double(jacobian(eq,vars));
    else
        warning('var_name not present in equation, proceeding without normalization')
    end
end
eq
disp('Coefficients:')
max_length = 1;
for j=1:length(vars)
    max_length = max(max_length,length(char(vars(j))));
end
for j=1:length(vars)
    txt = [' ',char(vars(j)),repmat(' ',1,2+max_length-length(char(vars(j)))+(sign(eq_jacob(j))==1)),...
            num2str(eq_jacob(j),'%.4f')];
    if roundn(eq_jacob(j),-4)==0
        txt = [txt,'  (',repmat(' ',1,(sign(eq_jacob(j))==1)),num2str(eq_jacob(j),'%.4e'),')'];
    end
    disp(txt)
end


