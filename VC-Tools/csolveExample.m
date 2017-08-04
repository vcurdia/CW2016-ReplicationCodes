% Example of a script to use csolve
%
% ..............................................................................
%
% Created: March 22, 2007 by Vasco Curdia
% Updated: July 26, 2011 by Vasco Curdia
% 
% Copyright 2007-2011 by Vasco Curdia

% -------------------------------------------------------------------------
%   Create function with system of equations
% -------------------------------------------------------------------------

% this function needs to accept as input a matrix x that will contain
% multiple candidates for the solution. the number of rows is the number of
% different variables being solved for, so that each column is a candidate.

% here I propose a way to create the function from within the script,
% instead o having to create "manually" the function. this is useful if the
% function is dependent on certain parameters from the original script.
% if all that the reader is looking for is an example of how the function
% looks like, just set the following option to 0 (zero).
DeleteFunction = 1;

% here the system of equations is very simple:
%   0 = y - a*log(z)
%   0 = z - b

a = 2;
b = 10;
vars = {'y','z'};
Nvars = length(vars);

fid=fopen('SysF.m','w');
fprintf(fid,'function f=SysF(x) \n');
fprintf(fid,'for j=1:size(x,2) \n'); % evaluate system of equations for each candidate
for jv=1:Nvars
    % for each variable (jv) in the system, fill it with the corresponding value 
    % in the current candidate (j)
    fprintf(fid,sprintf('    %s = x(%s,j); \n',vars{jv},int2str(jv)));
        % here the use of sprintf allows for the replacement of %s with a
        % string variable to be specified after the main argument
end
% write down each equation
fprintf(fid,sprintf('    f(1,j) = y - %f*log(z); \n',a));
fprintf(fid,sprintf('    f(2,j) = z - %f; \n',b));
    % the use of sprintf allows you to replace %f with a number,
    % specified after the main argument, in this latter case, b. The 
    % accuracy in this case is set to the default but if it needs to be
    % specified, use for example %.10 to specify precision at 1e-10.
fprintf(fid,'end  \n');
fclose(fid);

% -------------------------------------------------------------------------
%   Solve
% -------------------------------------------------------------------------

x0 = [1;1]; % the guess value
NumPrecision = 1e-10; % the precision of the solution for each equation

[xSol,rc] = csolve('SysF',x0,[],NumPrecision,1000);

% check to see if solution is not normal by looking at rc
if rc~=0, error(['Solution to system of equations is not normal, rc=',int2str(rc)]), end

% check the status of the system using the solution
disp(' ')
disp('System of equations if solution is plugged in:')
disp(SysF(xSol))

% present the results and set the variables to the solution:
disp(' ')
disp('Solution:')
for jv=1:Nvars
    disp(sprintf('%s = %f',vars{jv},xSol(jv)));
    eval(sprintf('%s = %f;',vars{jv},xSol(jv)))
end

% delete the function with the system of equations if it is not needed
% further
if DeleteFunction
    delete SysF.m
else
    open SysF
end
