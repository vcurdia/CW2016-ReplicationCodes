function pdflatex(fn,isCleanup)

% pdflatex
%
% Compiles tex file using PDFLaTeX.
%
% Usage:
%   pdflatex(fn)
%   pdflatex(fn,isCleanup)
%
% Inputs:
%
%   fn (string)
%   name of the tex file to compile, without extension
%
%   isCleanup (logical) [optional]
%   if set to 1, then all files with same file name but expensions other than
%   pdf or tex will be deleted. Default: 1
%
% ..............................................................................
%
% Created: August 22, 2011 by Vasco Curdia
% Updated: October 18, 2011 by Vasco Curdia
%          Does not show output on display. Can be retrieved from txt file
%          if isCleanup=0.
%
% Copyright 2011 by Vasco Curdia

%% -----------------------------------------------------------------------------

if ~exist('isCleanup','var'), isCleanup=1; end

eval(['!pdflatex ',fn,'.tex >> ',fn,'.txt'])
eval(['!pdflatex ',fn,'.tex >> ',fn,'.txt'])
eval(['!pdflatex ',fn,'.tex >> ',fn,'.txt'])

if isCleanup
  FileList = dir([fn,'.*']);
  FileList = {FileList(:).name};
  for j=1:length(FileList)
    if ~ismember(FileList{j}(end-2:end),{'tex','pdf'})
      delete(FileList{j})
    end
  end
end

%% -----------------------------------------------------------------------------
