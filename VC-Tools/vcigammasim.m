function [alpha,beta]=vcigammasim(ex,vx,bin)

% simulates the shape of a inverse-gamma distribution given the mean and
% variance
%
%   [alpha,beta]=vcigammasim(ex,vx,bin)
% 
% here I'm assuming the pdf of inverse-gamma to be:
% p(x;a,b)=b^(-a)*x^(-(a+1))*e^(-1/(b*x))/Gamma(a)
%
% ..............................................................................
%
% Created: March 10, 2004 by Vasco Curdia
% Updated: July 26, 2011 by Vasco Curdia
% 
% Copyright 2004-2011 by Vasco Curdia

%% ========================================================================

syms alpha beta
eq=[ex-1/beta/(alpha-1); 
    vx-1/(beta^2*(alpha-1)^2*(alpha-2))];

[alpha,beta]=solve(eq(1),eq(2),alpha,beta);
alpha=double(alpha);
beta=double(beta);
mode = 1/beta/(alpha+1);
simul=1./gamrnd(alpha,beta,100000,1);
figure
if nargin<3
    bin=100;
end
hist(simul,bin)
disp(['function: inv-gamma(',num2str(alpha),',',num2str(beta),')'])
disp(['mode: ',num2str(mode)])
bq = 1./gaminv([.95,.5,.05],alpha,beta);
disp(['5%, median and 95%: [ ',num2str(bq(1)),' , ',...
    num2str(bq(2)),' , ',num2str(bq(3)),' ]'])