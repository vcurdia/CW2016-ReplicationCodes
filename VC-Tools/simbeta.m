function simbeta(pmean,pse)

% Simulates the Beta distribution
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
a = pmean*(pmean-pmean^2-pse^2)/pse^2;
b = a*(1/pmean-1);
D.Dist = 'Beta';
D.Mean = pmean;
D.SE = pse;
D.a = a;
D.b = b;
D.Mode = (a-1)/(a+b-2);
for jprc=1:nprc
    eval(sprintf('D.prc%03.0f = betainv(%f,%.16f,%.16f);',...
        [1000,1]*Percentiles(jprc),a,b))
end
x = 0:1/nSample:1;
D.x = x;
D.pdf = betapdf(x,a,b);

%% ------------------------------------------------------------------------

%% Display results
disp(' ')
disp(D)
figure
plot(D.x,D.pdf,'-b')
xlim([0,1])
legend(D.Dist,'Location','NE')

%% ------------------------------------------------------------------------
