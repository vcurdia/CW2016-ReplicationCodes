function pdflatexall(varargin)

% pdflatexall
%
% Converts all tex files in directory to pdf using pdflatex
%
% .........................................................................
% 
% Created: January 13, 2010 by Vasco Curdia
% Updated: July 26, 2011 by Vasco Curdia
% Updated: August 22, 2011 by Vasco Curdia
%          Use pdflatex function.
% 
% Copyright 2010-2011 by Vasco Curdia

%% ------------------------------------------------------------------------

TEXFiles = dir('*.tex');
for j=1:length(TEXFiles)
    pdflatex(TEXFiles(j).name(1:end-4),varargin{:});
end

%% ------------------------------------------------------------------------
