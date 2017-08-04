function san=sadet(yy, freq,tsize)

% Conversion of GAUSS code for testing for seasonality
%
% GAUSS original from M. Watson
% Tests for seasonality
% Input:
%  y = Data series
%  freq= m - monthly, q - quarterly
%  tsize = size of test
% Output:
%  san = 0 no rejection
%        1 rejection
%        2 not enough obs
%
% ..............................................................................
%
% Created: October 29, 2002 by Vasco Curdia
% Updated: July 26, 2011 by Vasco Curdia
% 
% Copyright 2002-2011 by Vasco Curdia

if freq=='m'
    nf=12;
    rmin=48;
elseif freq=='q'
    nf=4;
    rmin=20;
else
    error('Data frequency must be either monthly (m) or quarterly (q)');
end

T=size(yy,1);
cvarsa=ones(T,1);
idm=fix(T/nf)+1;
sdm=repmat(eye(nf),idm,1);
sdm=[cvarsa sdm(1:size(yy,1),2:nf)]; % matrix of regressors
temp=[yy sdm];
temp(any(isnan(temp)'),:)=[];
if size(temp,1)<rmin
    san=2;
    return
end
yas=temp(:,1);
sdm=temp(:,2:size(temp,2));
xxi=rbinv(sdm'*sdm);
bsa=xxi*sdm'*yas;
err=yas-sdm*bsa;
s2hat=sum(err.^2)/(size(yas,1)-size(sdm,2));
vb=s2hat*xxi;
wstat=bsa(2:nf,:)'*inv(vb(2:nf,2:nf))*bsa(2:nf,:);
pv=1-chi2cdf(wstat,nf-1);
san=(pv<=tsize);
