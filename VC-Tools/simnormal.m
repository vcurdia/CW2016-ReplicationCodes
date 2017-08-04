function simnormal(pmean,pse)

% Simulates the Normal distribution
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
D.Dist = 'Normal';
D.Mean = pmean;
D.SE = pse;
D.Mode = pmean;
for jprc=1:nprc
    eval(sprintf('D.prc%03.0f = norminv(%f,%.16f,%.16f);',...
        [1000,1]*Percentiles(jprc),pmean,pse))
end
x = sort(normrnd(pmean,pse,nSample,1));
D.x = x;
D.pdf = normpdf(x,pmean,pse);

%% ------------------------------------------------------------------------

%% Display results
disp(' ')
disp(D)
figure
plot(D.x,D.pdf,'-b')
% xlim([D.prc010,D.prc990])
legend(D.Dist,'Location','NE')

%% ------------------------------------------------------------------------
