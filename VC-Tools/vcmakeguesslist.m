function []=vcmakeguesslist(sys,defaultguess)

% creates a text file with the list of variables in a system and presents
% them in such a way to create a structure with variable names and guess
% values that can then be updated and copied to another place
%
% required m-files:
%   - symbolic toolbox
%
% ..............................................................................
%
% Created: November 18, 2004 by Vasco Curdia
% Updated: July 26, 2011 by Vasco Curdia
% 
% Copyright 2004-2011 by Vasco Curdia

%--------------------------------------------------------------------------

fid=fopen('guesslist.txt','w');
vars = findsym(sys);
idx=findstr(',',vars);
nv = length(idx)+1;
for j=1:nv
    if j==1, vbegin=1; else, vbegin=idx(j-1)+2; end
    if j==nv, vend=length(vars); else, vend=idx(j)-1; end
    fprintf(fid,['guess(',num2str(j),').name  = ',vars(vbegin:vend),'; \n']);
    fprintf(fid,['guess(',num2str(j),').value = ',num2str(defaultguess),'; \n']);
end
fclose(fid);
open guesslist.txt