function D=siminvgamma(pmean,pse,isPlot)

% Compares the IG1 and IG2 distributions
%
% .........................................................................
% 
% Created: December 2, 2009 by Vasco Curdia
% Updated: July 26, 2011 by Vasco Curdia
% Updated: December 26, 2012 by Vasco Curdia
% 
% Copyright 2009-2012 by Vasco Curdia

%% ------------------------------------------------------------------------

%% Preamble
a0 = 5;
NumPrecision = 1e-10;
nMaxIt = 1000;
Percentiles = [0.01, 0.025, 0.05, 0.5, 0.95, 0.975, 0.99];
nprc = length(Percentiles);
nSample = 10000;
if ~exist('isPlot','var'),isPlot=1;end
%% ------------------------------------------------------------------------

%% IG1
if pse==inf
    a = 1;
else
    fname = sprintf('igamsolve%.0f',cputime*1e10);
    fid=fopen([fname,'.m'],'wt');
    fprintf(fid,'function f=%s(x)\n',fname);
    fprintf(fid,'pmean = %.16f;\n',pmean);
    fprintf(fid,'pvar = %.16f;\n',pse^2);
    fprintf(fid,'for j=1:length(x)\n');
    fprintf(fid,'    a = x(j);\n');
    fprintf(fid,'    f(j) = 1/(a-1)*(pmean*gamma(a)/gamma(a-1/2))^2-pmean^2-pvar;\n');
    fprintf(fid,'end\n');
    fclose(fid);
    [a,rc] = csolve(fname,a0,[],NumPrecision,nMaxIt);
    if rc~=0, error('Search for iGam parameters failed, rc=%.0f',rc), end
    delete([fname,'.m'])
end
b = (gamma(a-1/2)/pmean/gamma(a))^2;
D(1).Dist = 'IG1';
D(1).Mean = pmean;
D(1).SE = pse;
D(1).a = a;
D(1).b = b;
D(1).Mode = (1/b/(a+1/2))^(1/2);
for jprc=1:nprc
    eval(sprintf('D(1).prc%03.0f = gaminv(1-%f,%.16f,%.16f)^(-1/2);',...
        [1000,1]*Percentiles(jprc),a,b))
end
x = sort(gamrnd(a,b,nSample,1).^(-1/2));
D(1).x = x;
D(1).pdf = (x>0).*(gampdf(x.^(-2),a,b).*2./x.^3);

%% ------------------------------------------------------------------------

%% IG2
D(2).Dist = 'IG2';
if pse==inf
    a = 2;
else
    a = 2+pmean^2/pse^2;
end
b = 1/pmean/(a-1);
D(2).Mean = pmean;
D(2).SE = pse;
D(2).a = a;
D(2).b = b;
D(2).Mode = 1/b/(a+1);
for jprc=1:nprc
    eval(sprintf('D(2).prc%03.0f = gaminv(1-%f,%.16f,%.16f)^(-1);',...
        [1000,1]*Percentiles(jprc),a,b))
end
x = sort(gamrnd(a,b,nSample,1).^(-1));
D(2).x = x;
D(2).pdf = (x>0).*(gampdf(x.^(-1),a,b)./x.^2);

%% ------------------------------------------------------------------------

%% Display results
disp(' ')
disp(D(1))
disp(D(2))
if isPlot
  figure
  plot(D(1).x,D(1).pdf,'-b',D(2).x,D(2).pdf,'--r')
  xlim([0,max(D(1).prc990,D(2).prc990)])
  legend(D(1).Dist,D(2).Dist,'Location','NE')
end

%% ------------------------------------------------------------------------
