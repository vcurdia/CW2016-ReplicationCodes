function simgamma(pmean,pse)

% Simulates the Gamma distribution
%
% .........................................................................
% 
% Created: December 3, 2009 by Vasco Curdia
% Updated: July 26, 2011 by Vasco Curdia
% 
% Copyright 2009-2011 by Vasco Curdia

%% ------------------------------------------------------------------------

%% Preamble
a0 = 5;
NumPrecision = 1e-10;
nMaxIt = 1000;
Percentiles = [0.01, 0.025, 0.05, 0.5, 0.95, 0.975, 0.99];
nprc = length(Percentiles);
nSample = 10000;

%% ------------------------------------------------------------------------

%% Generate parameters and moments
a = (pmean/pse)^2;
b = pmean/a;
D.Dist = 'Gamma';
D.Mean = pmean;
D.SE = pse;
D.a = a;
D.b = b;
if a>=1
    D.Mode = (a-1)*b;
else
    D.Mode = NaN;
end
for jprc=1:nprc
    eval(sprintf('D.prc%03.0f = gaminv(%f,%.16f,%.16f);',...
        [1000,1]*Percentiles(jprc),a,b))
end
x = sort(gamrnd(a,b,nSample,1));
D.x = x;
D.pdf = gampdf(x,a,b);

%% ------------------------------------------------------------------------

%% Display results
disp(' ')
disp(D)
figure
plot(D.x,D.pdf,'-b')
% xlim([D.prc010,D.prc990])
legend(D.Dist,'Location','NE')

%% ------------------------------------------------------------------------
