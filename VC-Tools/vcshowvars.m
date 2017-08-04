function []=vcshowvars(vstruct,x,cb)

% vcshowvars
%
% Presents each variable of the structure according to its name in the 
% structure and evaluates its value, comparing to the benchmark
%
% []=vcshowvars(vstruct,x,cb)
%
% inputs:
%  vstruct - a structure with at least two fields:
%       name - containing the name of the variables to show
%       benchmark - containing the benchmark number for the variables
%  x - values of variables
%  cb - binary variable to denote whether to compare to benchmark
%
% .........................................................................
%
% Created: April 09, 2005 by Vasco Curdia
% Updated: July 26, 2011 by Vasco Curdia
% 
% Copyright 2005-2011 by Vasco Curdia

% -------------------------------------------------------------------------

if nargin == 2
    cb = 0;
end
nV = length(vstruct);
maxnmlength = length(vstruct(1).name);
for j=2:nV
    maxnmlength = max(maxnmlength,length(vstruct(j).name));
end
for j=1:nV
    nmtext = vstruct(j).name;
    ldiff = maxnmlength-length(nmtext);
    if ldiff>0, nmtext = [nmtext, repmat(' ',1,ldiff)]; end
    if sign(x(j))~=-1, vnntxt = ' '; else vnntxt = ''; end
    if cb
        if vstruct(j).benchmark == 0
            disp(sprintf([nmtext,' = ',vnntxt,'%f (------)'],x(j)))
        else
            if sign(x(j)/vstruct(j).benchmark)~=-1, cbnntxt = ' '; else cbnntxt = ''; end
            disp(sprintf([nmtext,' = ',vnntxt,'%f (',cbnntxt,'%.3f)'],x(j),x(j)/vstruct(j).benchmark))
        end
    else
        disp(sprintf([nmtext,' = ',vnntxt,'%f'],x(j)))
    end
end
