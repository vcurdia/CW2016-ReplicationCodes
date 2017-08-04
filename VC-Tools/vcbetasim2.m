function [alpha,betta]=vcbetasim2(mx,vx,bin)

% Simulates the shape of a beta distribution given the mode and variance
%
%   [alpha,betta]=vcbetasim2(mx,vx,bin)
%
% ..............................................................................
%
% Created by Vasco Curdia in March 2004
% Modified by Daria Finocchiaro in August 2004
% Updated: July 26, 2011 by Vasco Curdia
% 
% Copyright 2004-2011 by Vasco Curdia and Daria Finocchiaro

%% ========================================================================

syms alpha betta NA
eq=[mx-((alpha-1)/(alpha+betta-2)); 
    vx-alpha*betta/((alpha+betta)^2*(alpha+betta+1))];
[alpha,betta]=solve(eq(1),eq(2),alpha,betta);

alpha=double(alpha);
betta=double(betta);
n = size(alpha,1);

% If beta and alpha are both greater than 1 the beta distribution has a
% unique mode given by mx=(alpha-1)/(alpha+beta-2)
for j=1:n
    a=alpha(j,:);
    b=betta(j,:);
    if imag(a)<1.e-016, a=real(a); end
    if imag(b)<1.e-016, b=real(b); end
    if [a b]>[1 1]
        ex=a/(a+b);
        a;
        b;
        simul=betarnd(a,b,1e+5,1);
        figure
        if nargin<3
            bin=100;
        end
        hist(simul,bin)
        disp(['function: beta(',num2str(a),',',num2str(b),')'])
        disp(['mean: ',num2str(ex)])
        bq = betainv([.025,.5,.975],a,b);
        disp(['2.5%, median and 97,5%: [ ',num2str(bq(1)),' , ',...
                num2str(bq(2)),' , ',num2str(bq(3)),' ]'])
        break
    else disp('Not Unimodal')
        ex=NA
    end
end
alpha = a;
betta = b;