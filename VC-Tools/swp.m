function swp(A,varargin)

% SWP
%
% converts a matrix to a table in SWP
%
% Usage:
%   swp(A)
%   swp(A,FileName)
%   swp(A,NumPrecision)
%   swp(A,FileName,NumPrecision)
%   swp(A,NumPrecision,FileName)
%
% Inputs:
%   A [double]
%       matrix variable already existing in the workspace
%   filename [string] [Optional]
%       name of the filename, excluding ".tex" suffix
%       Default: 'temp'
%   NumDec [double] [Optional]
%       Fixed decimal. Default: 4 
%
% Description:
% This function creates a new file which is readable in both latex and
% SWP. This file is created in the current working directory.
%
% .........................................................................
%
% Created: March 11, 2003 by Vasco Curdia
% 
% Copyright 2003-2017 by Vasco Curdia

%% ------------------------------------------------------------------------

%% check arguments
if nargin==1
    FileName = 'temp';
    NumDec = 4;
elseif nargin==2
    if ischar(varargin{1})
        FileName = varargin{1};
        NumDec = 4;
    else
        NumDec = varargin{1};
        FileName = 'temp';
    end
elseif nargin==3
    if ischar(varargin{1})
        FileName = varargin{1};
        NumDec = varargin{2};
    else
        FileName = varargin{2};
        NumDec = varargin{1};
    end
else
    error('Too many arguments')
end

%% Create tex file
[nr,nc]=size(A);
cc=repmat('c',1,nc);
fid=fopen([FileName,'.tex'],'w');
fprintf(fid,'\n\\documentclass{article}\n');
% fprintf(fid,'%%TCIDATA{OutputFilter=LATEX.DLL} \n');
% fprintf(fid,'%%TCIDATA{Version=4.00.0.2321} \n');
% fprintf(fid,'%%TCIDATA{Created=Tuesday, March 11, 2003 22:13:46} \n');
% fprintf(fid,'%%TCIDATA{LastRevised=Tuesday, March 11, 2003 22:32:46} \n');
% fprintf(fid,'%%TCIDATA{<META NAME="GraphicsSave" CONTENT="32">} \n');
% fprintf(fid,'%%TCIDATA{<META NAME="DocumentShell" CONTENT="Standard LaTeX\\Blank - Standard LaTeX Article">} \n');
% fprintf(fid,'%%TCIDATA{CSTFile=40 LaTeX article.cst} \n');
% fprintf(fid,'\\input{tcilatex}\n');
fprintf(fid,'\\begin{document}\n');
fprintf(fid,'\\begin{tabular}{%s}\n',cc);
for i=1:nr
    text=sprintf(['%.',int2str(NumDec),'f'],A(i,1));
    for j=2:nc
        text=sprintf(['%s & %.',int2str(NumDec),'f'],text,A(i,j));
    end
    if i==nr
        text=[text,'\n'];
    else
        text=[text,'\\\\\n'];
    end
    fprintf(fid,text);
end
fprintf(fid,'\\end{tabular}\n');
fprintf(fid,'\\end{document}\n');
fclose(fid);

%% ------------------------------------------------------------------------
