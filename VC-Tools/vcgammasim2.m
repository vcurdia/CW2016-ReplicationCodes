function [alpha,beta]=vcgammasim2(mx,vx,bin)

% simulates the shape of a gamma distribution given the mode and variance
%
%   [alpha,beta]=vcgammasim2(mx,vx,bin)
%
% Pdf of Gamma :
% p(x;a,b)=b^(-a)*x^(a-1)*e^(-x/b)/Gamma(a)
%
% ..............................................................................
%
% Created: September 2004 by Daria Finocchiaro
% Updated: July 26, 2011 by Vasco Curdia
% 
% Copyright 2004-2011 by Vasco Curdia and Daria Finocchiaro

%--------------------------------------------------------------------------

syms alpha beta
eq=[mx-beta*(alpha-1); 
    vx-alpha*beta^2];
[alpha,beta]=solve(eq(1),eq(2),alpha,beta);
alpha=double(alpha);
beta=double(beta);
n = size(alpha,1);
for j=1:n
    a=alpha(j,:);
    b=beta(j,:);
    if imag(a)<1.e-030, a=real(a); end
    if imag(b)<1.e-030, b=real(b); end
    if a>=1
        ex=b*a;
        a;
        b;
        simul=gamrnd(a,b,1e+5,1);
        figure
        if nargin<3
            bin=100;
        end
        hist(simul,bin)
        disp(['function: gamma(',num2str(a),',',num2str(b),')'])
        disp(['mean: ',num2str(ex)])
        bq = gaminv([.025,.5,.975],a,b);
        disp(['2.5%, median and 97,5%: [ ',num2str(bq(1)),' , ',...
                num2str(bq(2)),' , ',num2str(bq(3)),' ]'])
        break
    end
end
alpha = a;
beta = b;