function [alpha,betta]=vcgammasim(ex,vx,bin)

% simulates the shape of a gamma distribution given the mean and variance
%
%   [alpha,betta]=vcgammasim(ex,vx,bin)
%
% Pdf of Gamma :
% p(x;a,b)=b^(-a)*x^(a-1)*e^(-x/b)/Gamma(a)
%
% ..............................................................................
%
% Created: January 27, 2005 by Vasco Curdia
% Updated: July 26, 2011 by Vasco Curdia
% 
% Copyright 2005-2011 by Vasco Curdia

%% ========================================================================

syms alpha betta
eq=[ex-betta*alpha; 
    vx-alpha*betta^2];
[alpha,betta]=solve(eq(1),eq(2),alpha,betta);
alpha=double(alpha);
betta=double(betta);
n = size(alpha,1);
for j=1:n
    a=alpha(j,:);
    b=betta(j,:);
    if imag(a)<1.e-030, a=real(a); end
    if imag(b)<1.e-030, b=real(b); end
    if a>=1
        mx=b*(a-1);
        a;
        b;
        simul=gamrnd(a,b,1e+5,1);
        figure
        if nargin<3
            bin=100;
        end
        hist(simul,bin)
        disp(['function: gamma(',num2str(a),',',num2str(b),')'])
        disp(['mode: ',num2str(mx)])
        bq = gaminv([.025,.5,.975],a,b);
        disp(['2.5%, median and 97,5%: [ ',num2str(bq(1)),' , ',...
                num2str(bq(2)),' , ',num2str(bq(3)),' ]'])
        break
    end
end
alpha = a;
betta = b;