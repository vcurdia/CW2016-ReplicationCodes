function fid = createtex(fName,fTitle)

% createtex
%
% .........................................................................
%
% Created: February, 2016 by Vasco Curdia
% 
% Copyright 2016-2017 by Vasco Curdia

% -------------------------------------------------------------------------

fid = fopen(sprintf('%s.tex',fName),'wt');

fprintf(fid,'\\documentclass[12pt]{article}\n\n');

if ~exist('fTitle','var'), fTitle = ''; end
fTitle = strrep(fTitle,'_',' ');

fprintf(fid,'\\newcommand{\\MyTitle}{%s}\n',fTitle);
fprintf(fid,'\\newcommand{\\MyShortTitle}{%s}\n\n',strrep(fTitle,'\\',' '));

fprintf(fid,'\\usepackage{cmap}\n');
fprintf(fid,'\\usepackage{amsfonts}\n');
fprintf(fid,'\\usepackage{amssymb}\n');
fprintf(fid,'\\usepackage{amsmath}\n');
fprintf(fid,'\\usepackage{graphicx}\n');
fprintf(fid,'\\usepackage{color}\n\n');

fprintf(fid,'\\usepackage[\n');
fprintf(fid,'  pdftex,\n');
fprintf(fid,'  pdfstartview=Fit,\n');
fprintf(fid,'  pdfpagelayout=SinglePage,\n');
fprintf(fid,'  bookmarksopen,\n');
fprintf(fid,'  bookmarksopenlevel=1,\n');
fprintf(fid,'  colorlinks,\n');
fprintf(fid,'  linkcolor=RefColor,\n');
fprintf(fid,'  citecolor=RefColor,\n');
fprintf(fid,'  urlcolor=RefColor,\n');
fprintf(fid,'  filecolor=RefColor,\n');
fprintf(fid,'  naturalnames\n');
fprintf(fid,']{hyperref}\n');
fprintf(fid,...
        '\\hypersetup{pdftitle={\\MyShortTitle}}\n\n');

%fprintf(fid,'\\usepackage{authblk}\n\n');

fprintf(fid,'\\usepackage[letterpaper]{geometry}\n');
fprintf(fid,'\\geometry{top=1.25in, bottom=1.0in, left=1.0in, right=1.0in}\n');
fprintf(fid,'\\setlength{\\headheight}{15pt}\n\n');

fprintf(fid,'\\usepackage{setspace}\n\n');

fprintf(fid,'\\usepackage{fancyhdr}\n');
fprintf(fid,'\\pagestyle{fancy}\n');
fprintf(fid,'\\fancyhf{}\n');
fprintf(fid,'\\fancyhead[C]{\\textsc{\\MyShortTitle}}\n');
fprintf(fid,'\\fancyfoot[C]{\\thepage}\n');
fprintf(fid,'\\renewcommand{\\headrulewidth}{0pt}\n\n');

fprintf(fid,'\\definecolor{RefColor}{rgb}{0.1,0.3,0.5}\n\n');

fprintf(fid,'\\begin{document}\n\n');

fprintf(fid,'\\title{\\MyTitle}\n\n');

fprintf(fid,'\\date{\\today}\n\n');

fprintf(fid,'\\maketitle\n\n');

fprintf(fid,'\\thispagestyle{empty}\n\n');

