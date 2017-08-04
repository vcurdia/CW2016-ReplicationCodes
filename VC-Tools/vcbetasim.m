function [alpha,betta]=vcbetasim(ex,vx,bin)
% simulates the shape of a beta distribution given the mean and variance
%
% ..............................................................................
%
% Created: March 06, 2004 by Vasco Curdia
% Updated: July 26, 2011 by Vasco Curdia
% 
% Copyright 2004-2011 by Vasco Curdia

%--------------------------------------------------------------------------

syms alpha betta
eq=[ex-alpha/(alpha+betta); 
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
        bmode = (a-1)/(a+b-2);
        a;
        b;
        simul=betarnd(a,b,1e+5,1);
        figure
        if nargin<3
            bin=100;
        end
        hist(simul,bin)
        disp(['function: beta(',num2str(a),',',num2str(b),')'])
        disp(['mode: ',num2str(bmode)])
        bq = betainv([.025,.5,.975],a,b);
        disp(['2.5%, median and 97,5%: [ ',num2str(bq(1)),' , ',...
                num2str(bq(2)),' , ',num2str(bq(3)),' ]'])
        break
    else 
        disp('Not Unimodal')
    end
end
alpha = a;
betta = b;