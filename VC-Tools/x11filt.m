function x11=x11filt(freq)

% Conversion of GAUSS code for computing the X11 filter
% GAUSS original from M. Watson
% x11filt.prc, mww 3/29/00
% Compute X11 Filter
% Follow 8 Steps in Watson's JBES discussion of Ghysels, Granger, Siklos
% (Which is taken from Wallis's 1974 JASA paper) 
%
% ..............................................................................
%
% Created: October 29, 2002 by Vasco Curdia
% Updated: July 26, 2011 by Vasco Curdia
% 
% Copyright 2002-2011 by Vasco Curdia

if freq=='m'
    nf=12;
elseif freq=='q'
    nf=4;
else
    error('Data frequency must be either monthly (m) or quarterly (q)');
end


% Step 1: TC1 = a1(L)x(t)
  a1=zeros(nf+1,1);
  a1(1)=1/(2*nf);
  a1(nf+1)=1/(2*nf);
  a1(2:nf)=1/nf*ones(nf-1,1);

% Step 2: SI1=x-TC1=a2(L)x 
  a2=-a1;
  a2(nf/2+1)=1+a2(nf/2+1);

% Step 3: S1=S1(L)SI1=a3(L)x 
  s1=zeros(4*nf+1,1);
  s1(1)=1/9;
  s1(nf+1)=2/9;
  s1(2*nf+1)=3/9;
  s1(3*nf+1)=2/9;
  s1(4*nf+1)=1/9;

  temp=size(s1,1)-size(a2,1);
  temp=temp/2;
  temp=zeros(temp,1);
  temp=[temp;a2;temp];

  a3=conv(s1,temp);
  a3=a3(find(a3));

% Step 4: S2=a2(L)S1=a4(L)x 
  temp=size(a3,1)-size(a2,1);
  temp=temp/2;
  temp=zeros(temp,1);
  temp=[temp;a2;temp];
  a4=conv(a3,temp);
  a4=a4(find(a4));

% Step 5: TC2=H(L)(x-S2)=a5(L)x 
  h=zeros(nf+1,1);
  if freq=='m'
      h(1)=-.0194;
      h(2)=-.0279;
      h(3)=0;
      h(4)=.0655;
      h(5)=.1474;
      h(6)=.2143;
      h(7)=.2402;
      h(8:13)=flipud(h(1:6));
  elseif freq=='q'
      h(1)=-.073;
      h(2)=.294;
      h(3)=.558;
      h(4:5)=flipud(h(1:2));
  end

  temp1=-a4;
  mid=size(temp1,1)-1;
  mid=mid/2;
  mid=mid+1;
  temp1(mid)=1+temp1(mid);

  temp=size(temp1,1)-size(h,1);
  temp=temp/2;
  temp=zeros(temp,1);
  temp=[temp;h;temp];
  a5=conv(temp1,temp);
  a5=a5(find(a5));

% Step 6: S3=S3(L)(x-TC2)=a6(L)x 
  temp1=-a5;
  mid=size(temp1,1)-1;
  mid=mid/2;
  mid=mid+1;
  temp1(mid)=1+temp1(mid);
  s3=zeros(6*nf,1);
  s3(1)=1/15;
  s3(nf+1)=2/15;
  s3(2*nf+1)=3/15;
  s3(3*nf+1)=3/15;
  s3(4*nf+1)=3/15;
  s3(5*nf+1)=2/15;
  s3(6*nf+1)=1/15;
  temp=size(temp1,1)-size(s3,1);
  temp=temp/2;
  temp=zeros(temp,1);
  temp=[temp;s3;temp];
  a6=conv(temp1,temp);
  a6=a6(find(a6));

% Step 7: S4=a2(L)S3=a7(L)x 
  temp=size(a6,1)-size(a2,1);
  temp=temp/2;
  temp=zeros(temp,1);
  temp=[temp;a2;temp];
  a7=conv(a6,temp);
  a7=a7(find(a7));

% Step 8: XSA=x-S4=a8x 
  a8=-a7;
  mid=size(a8,1)-1;
  mid=mid/2;
  mid=mid+1;
  a8(mid)=1+a8(mid);
  x11=a8(find(a8));
