% testrecessionshades
%
% Test the function to put recession shades.
%
% ..............................................................................
%
% Created: April 26, 2011 by Vasco Curdia
% Updated: July 26, 2011 by Vasco Curdia
% 
% Copyright 2011 by Vasco Curdia


%% Prepare information
x = 1:10;
y = -2+0.5*(1:10);
tid = 4*(1:10);
Rec = [zeros(1,2),ones(1,2),zeros(1,3),ones(1,2),0];

%% Make plots
plot(tid,x,'r')
hold on
plot(tid,y,'b')
hold off
legend('x','y','Location','NW')

%% Add Recession Shades
recessionshades(Rec,'TimeIdx',tid)

%% Make PDF
% If want to test how it looks in pdf, and have epstopdf installed
% print -depsc testplot.eps
% !epstopdf testplot.eps
% winopen testplot.pdf
