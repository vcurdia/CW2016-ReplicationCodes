function epstopdfall

% epstopdfall
%
% Converts all eps files in directory to pdf using epstopdf
%
% .........................................................................
% 
% Created: July 21, 2009 by Vasco Curdia
% Updated: July 26, 2011 by Vasco Curdia
% 
% Copyright 2009-2011 by Vasco Curdia

%% ------------------------------------------------------------------------

EPSFiles = dir('*.eps');
for j=1:length(EPSFiles)
    system(['epstopdf ',EPSFiles(j).name]);
end

%% ------------------------------------------------------------------------
