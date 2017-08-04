function [G1,C,impact,fmat,fwt,ywt,gev,eu]=fastgensysJaeWonvb(...
  g0,g1,c,psi,pi,fid,verbose,div,realsmall,UsePinv)

% fastgensysJaeWonvb
%
% based on fastgensysJaeWon, but added options for display manipulation
%
% Options added:
%   fid
%     If 1, then output shown to screen. Otherwise it is the handle of a text
%     file currently open. Default: 1
%   verbose (logical)
%     If set to 1 the displays certain messages. Default: 1.
%   div
%     value of div to be used, instead of generated.
%     already present in original gensys, but it's now mixed w/ other args
%   realsmall
%     can choose the level of realsmall. default 1e-8 (CS: 1e-6)
%   UsePinv
%     whether to use pinv or inv to invert G0. 
%     in some cases it is better to make realsmall smaller and use instead pinv
%     default: 1 (CS: 0)
%
% ------------------------------------------------------------------------------
% function [G1,C,impact,fmat,fwt,ywt,gev,eu]=gensys(g0,g1,c,psi,pi,div)
% System given as
%        g0*y(t)=g1*y(t-1)+c+psi*z(t)+pi*eta(t),
% with z an exogenous variable process and eta being endogenously determined
% one-step-ahead expectational errors.  Returned system is
%       y(t)=G1*y(t-1)+C+impact*z(t)+ywt*inv(I-fmat*inv(L))*fwt*z(t+1) .
% If z(t) is i.i.d., the last term drops out.
% If div is omitted from argument list, a div>1 is calculated.
% eu(1)=1 for existence, eu(2)=1 for uniqueness.  eu(1)=-1 for
% existence only with not-s.c. z; eu=[-2,-2] for coincident zeros.
% By Christopher A. Sims
% Corrected 10/28/96 by CAS
% ------------------------------------------------------------------------------
% 
% Created: July 28, 2011 by Vasco Curdia
% Updated: July 14, 2014 by Vasco Curdia
% 
% Copyright 2011-2014 by Vasco Curdia

if nargin<6 || isempty(fid), fid = 1; end
if nargin<7 || isempty(verbose), verbose = 1; end

fixdiv = exist('div','var') && ~isempty(div);

if ~exist('realsmall','var') || isempty(realsmall)
    realsmall = 1e-8;
end

if ~exist('UsePinv','var') || isempty(UsePinv)
    UsePinv = 1;
end

eu=[0;0];
n=size(g0,1);
[a b q z v]=qz(g0,g1);
if ~fixdiv, div=1.01; end
nunstab=0;
zxz=0;
for i=1:n
% ------------------div calc------------
   if ~fixdiv
      if abs(a(i,i)) > 0
         divhat=abs(b(i,i))/abs(a(i,i));
	 % bug detected by Vasco Curdia and Daria Finocchiaro, 2/25/2004  A root of
	 % exactly 1.01 and no root between 1 and 1.02, led to div being stuck at 1.01
	 % and the 1.01 root being misclassified as stable.  Changing < to <= below fixes this.
         if 1+realsmall<divhat && divhat<=div
            div=.5*(1+divhat);
         end
      end
   end
% ----------------------------------------
   nunstab=nunstab+(abs(b(i,i))>div*abs(a(i,i)));
   if abs(a(i,i))<realsmall && abs(b(i,i))<realsmall
      zxz=1;
   end
end
% div;
% 
% nunstab;
if ~zxz
   [a b q z]=qzdiv_ordqz(div,a,b,q,z);  %[a b q z]=qzdiv(div,a,b,q,z);   
end
gev=[diag(a) diag(b)];
if zxz
   if verbose
     fprintf(fid,'Coincident zeros. Indeterminacy and/or nonexistence.\n');
   end
   eu=[-2;-2];
   % correction added 7/29/2003.  Otherwise the failure to set output
   % arguments leads to an error message and no output (including eu).
   G1=[];C=[];impact=[];fmat=[];fwt=[];ywt=[];gev=[];
   return
end
q1=q(1:n-nunstab,:);
q2=q(n-nunstab+1:n,:);
z1=z(:,1:n-nunstab)';
z2=z(:,n-nunstab+1:n)';
a2=a(n-nunstab+1:n,n-nunstab+1:n);
b2=b(n-nunstab+1:n,n-nunstab+1:n);
etawt=q2*pi;
% zwt=q2*psi;
[ueta,deta,veta]=svd(etawt);
md=min(size(deta));
bigev=find(diag(deta(1:md,1:md))>realsmall);
ueta=ueta(:,bigev);
veta=veta(:,bigev);
deta=deta(bigev,bigev);
% ------ corrected code, 3/10/04
%bigev
%nunstab
%pause

eu(1) = length(bigev)>=nunstab;
if ~eu(1)==1 && verbose
    fprintf(fid,'eu(1): %.0f, bigev: %.0f, nunstab: %.0f, npi: %.0f\n',...
            eu(1),length(bigev),nunstab,size(pi,2));
end    
% ------ Code below allowed "existence" in cases where the initial lagged state was free to take on values
% ------ inconsistent with existence, so long as the state could w.p.1 remain consistent with a stable solution
% ------ if its initial lagged value was consistent with a stable solution.  This is a mistake, though perhaps there
% ------ are situations where we would like to know that this "existence for restricted initial state" situation holds.
%% [uz,dz,vz]=svd(zwt);
%% md=min(size(dz));
%% bigev=find(diag(dz(1:md,1:md))>realsmall);
%% uz=uz(:,bigev);
%% vz=vz(:,bigev);
%% dz=dz(bigev,bigev);
%% if isempty(bigev)
%% 	exist=1;
%% else
%% 	exist=norm(uz-ueta*ueta'*uz) < realsmall*n;
%% end
%% if ~isempty(bigev)
%% 	zwtx0=b2\zwt;
%% 	zwtx=zwtx0;
%% 	M=b2\a2;
%% 	for i=2:nunstab
%% 		zwtx=[M*zwtx zwtx0];
%% 	end
%% 	zwtx=b2*zwtx;
%% 	[ux,dx,vx]=svd(zwtx);
%% 	md=min(size(dx));
%% 	bigev=find(diag(dx(1:md,1:md))>realsmall);
%% 	ux=ux(:,bigev);
%% 	vx=vx(:,bigev);
%% 	dx=dx(bigev,bigev);
%% 	existx=norm(ux-ueta*ueta'*ux) < realsmall*n;
%% else
%% 	existx=1;
%% end
% ----------------------------------------------------
% Note that existence and uniqueness are not just matters of comparing
% numbers of roots and numbers of endogenous errors.  These counts are
% reported below because usually they point to the source of the problem.
% ------------------------------------------------------
[ueta1,deta1,veta1]=svd(q1*pi);
md=min(size(deta1));
bigev=find(diag(deta1(1:md,1:md))>realsmall);
ueta1=ueta1(:,bigev);
veta1=veta1(:,bigev);
deta1=deta1(bigev,bigev);
%% if existx | nunstab==0
%%    %disp('solution exists');
%%    eu(1)=1;
%% else
%%     if exist
%%         %disp('solution exists for unforecastable z only');
%%         eu(1)=-1;
%%     %else
%%         %fprintf(1,'No solution.  %d unstable roots. %d endog errors.\n',nunstab,size(ueta1,2));
%%     end
%%     %disp('Generalized eigenvalues')
%%    %disp(gev);
%%    %md=abs(diag(a))>realsmall;
%%    %ev=diag(md.*diag(a)+(1-md).*diag(b))\ev;
%%    %disp(ev)
%% %   return;
%% end
if isempty(veta1)
	unique=1;
else
	unique=norm(veta1-veta*veta'*veta1)<realsmall*n;
end
if unique
   %disp('solution unique');
   eu(2)=1;
else
   if verbose
     fprintf(fid,'Indeterminacy. %d loose endog errors.\n',size(veta1,2)-size(veta,2));
   end
   %disp('Generalized eigenvalues')
   %disp(gev);
   %md=abs(diag(a))>realsmall;
   %ev=diag(md.*diag(a)+(1-md).*diag(b))\ev;
   %disp(ev)
%   return;
end
tmat = [eye(n-nunstab) -(ueta*(deta\veta')*veta1*deta1*ueta1')'];
G0= [tmat*a; zeros(nunstab,n-nunstab) eye(nunstab)];
G1= [tmat*b; zeros(nunstab,n)];
% ----------------------
% G0 is always non-singular because by construction there are no zeros on
% the diagonal of a(1:n-nunstab,1:n-nunstab), which forms G0's ul corner.
% -----------------------
if UsePinv
    G0I = pinv(G0);
else
    G0I = inv(G0);
end
G1=G0I*G1;
usix=n-nunstab+1:n;
C=G0I*[tmat*q*c;(a(usix,usix)-b(usix,usix))\q2*c];
impact=G0I*[tmat*q*psi;zeros(nunstab,size(psi,2))];
fmat=b(usix,usix)\a(usix,usix);
fwt=-b(usix,usix)\q2*psi;
ywt=G0I(:,usix);
% -------------------- above are output for system in terms of z'y -------
G1=real(z*G1*z');
C=real(z*C);
impact=real(z*impact);
% Correction 10/28/96:  formerly line below had real(z*ywt) on rhs, an error.
ywt=z*ywt;
