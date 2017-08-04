function ypad=padar(y,n,arvec)

% Conversion of GAUSS code to pad data series
%
% GAUSS original from M. Watson
% Pad Data series y out using AR Forecasts and Backcasts
%     y -- series to be padded
%     n -- number of terms to pad forward and backward
%     arvec -- vector of AR lags
%             if lags 1,3 and 6 are needed, then ARVEC=1|3|6, etc.
%
% ..............................................................................
%
% Created: October 29, 2002 by Vasco Curdia
% Updated: July 26, 2011 by Vasco Curdia
% 
% Copyright 2002-2011 by Vasco Curdia


nar=max(arvec);

% Pad out future
w=y(nar+1:size(y,1));
x=ones(size(w,1),1);
for i=1:size(arvec,1)
    x=[x y(nar+1-arvec(i):size(y,1)-arvec(i))];
end
bols=inv(x'*x)*(x'*w);
beta=zeros(1+nar,1);
beta(1)=bols(1);
for i=1:size(arvec,1)
    beta(1+arvec(i))=bols(i+1);
end
v=flipud(y(size(y,1)-nar+1:size(y,1)));
forc=zeros(n,1);
for i=1:n
    forc(i)=beta'*[1;v];
    v(2:size(v,1))=v(1:size(v,1)-1);
    v(1)=forc(i);
end
ypad=[y;forc];

% Pad out past, by reversing series 
y=flipud(y);
w=y(nar+1:size(y,1));
x=ones(size(w,1),1);
for i=1:size(arvec,1)
    x=[x y(nar+1-arvec(i):size(y,1)-arvec(i))];
end
bols=inv(x'*x)*(x'*w);
beta=zeros(1+nar,1);
beta(1)=bols(1);
for i=1:size(arvec,1)
    beta(1+arvec(i))=bols(i+1);
end
v=flipud(y(size(y,1)-nar+1:size(y,1)));
forc=zeros(n,1);
for i=1:n
    forc(i)=beta'*[1;v];
    v(2:size(v,1))=v(1:size(v,1)-1);
    v(1)=forc(i);
end
forc=flipud(forc);
ypad=[forc;ypad];
