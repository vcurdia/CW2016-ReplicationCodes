function [alpha,beta]=vcigammasim2(mx,vx,bin)

% simulates the shape of a gamma distribution given the mode and variance
%
%   [alpha,beta]=vcigammasim2(mx,vx,bin)
%
% here I'm assuming the pdf of i-gamma to be:
% p(x;a,b)=b^(-a)*x^(-(a+1))*e^(-1/(b*x))/Gamma(a)
%
% ..............................................................................
%
% Created: March 10, 2004 by Vasco Curdia
% Updated: September 2004 by Daria Finocchiaro
% Updated: July 26, 2011 by Vasco Curdia
% 
% Copyright 2004-2011 by Vasco Curdia

%% ========================================================================

syms alpha beta
eq=[mx-1/beta/(alpha+1); 
    vx-1/(beta^2*(alpha-1)^2*(alpha-2))];
[alpha,beta]=solve(eq(1),eq(2),alpha,beta);
alpha=double(alpha);
beta=double(beta);
n = size(alpha,1);
for j=1:n
    a=alpha(j,:);
    b=beta(j,:);
    if imag(a)<1.e-030, a=real(a); end
    if imag(b)<1.e-030, b=real(b); end
    if a>2
        ex=1/b/(a-1);
        a;
        b;
        simul=1./gamrnd(a,b,1e+5,1);
        figure
        if nargin<3
            bin=200;
        end
        hist(simul,bin)
        disp(['function: inv-gamma(',num2str(a),',',num2str(b),')'])
        disp(['mean: ',num2str(ex)])
        bq = 1./gaminv([.975,.5,.025],a,b);
        disp(['2.5%, median and 97,5%: [ ',num2str(bq(1)),' , ',...
                num2str(bq(2)),' , ',num2str(bq(3)),' ]'])
        break
    end
end
alpha = a;
beta = b;