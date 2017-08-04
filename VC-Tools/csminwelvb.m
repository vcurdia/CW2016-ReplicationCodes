function [fh,xh,gh,H,itct,fcount,retcodeh] = csminwelvb(fcn,x0,H0,grad,crit,nit,verbose,fn,varargin)

% csminwelvb
%
% builds on csminwel from Christopher A. Sims but makes a few tweaks
%
% Original Description
% %[fhat,xhat,ghat,Hhat,itct,fcount,retcodehat] = csminwel(fcn,x0,H0,grad,crit,nit,varargin)
% % fcn:   string naming the objective function to be minimized
% % x0:    initial value of the parameter vector
% % H0:    initial value for the inverse Hessian.  Must be positive definite.
% % grad:  Either a string naming a function that calculates the gradient, or the null matrix.
% %        If it's null, the program calculates a numerical gradient.  In this case fcn must
% %        be written so that it can take a matrix argument and produce a row vector of values.
% % crit:  Convergence criterion.  Iteration will cease when it proves impossible to improve the
% %        function value by more than crit.
% % nit:   Maximum number of iterations.
% % varargin: A list of optional length of additional parameters that get handed off to fcn each
% %        time it is called.
% %        Note that if the program ends abnormally, it is possible to retrieve the current x,
% %        f, and H from the files g1.mat and H.mat that are written at each iteration and at each
% %        hessian update, respectively.  (When the routine hits certain kinds of difficulty, it
% %        write g2.mat and g3.mat as well.  If all were written at about the same time, any of them
% %        may be a decent starting point.  One can also start from the one with best function value.)
%
% Added option: verbose (double)
% Allows control of the amount of information displayed:
%  0: no output at all displayed
%  1: only minimal information shown: it count, f level and improvement
%  2: full display as in original csminwel
%
% Added additional input/output argument: fn
% if fn is a nonempty string then verbose is routed to the filename with name fn.
% fn can also be the handle for a file opened in the calling program.
% Otherwise verbose is shown on screen.  Default: []
%
% ..............................................................................
%
% Created: February 14, 2011 by Vasco Curdia
%          - Relative to csminwel converted logical operators into
%            short-circuited logical operators.
%          - do not save the g1, g2, g3 matrices.
%          - added Verbose to list of options in the function.
% Updated: July 26, 2011 by Vasco Curdia
% Updated: July 28, 2011 by Vasco Curdia
%          changed the way things are printed
% Updated: July 2, 2012 by Vasco Curdia
%          uses now numgradvb
%
% Copyright 2011-2012 by Vasco Curdia

[nx,no]=size(x0);
nx=max(nx,no);
if ~exist('verbose','var') || isempty(verbose), verbose = 1; end
if ~exist('fn','var') || isempty(fn)
  fid=1;
elseif ischar(fn)
  fid=fopen(fn,'wt');
else
  fid = fn;
end
NumGrad = isempty(grad);
done=0;
itct=0;
fcount=0;
snit=100;
%tailstr = ')';
%stailstr = [];
% Lines below make the number of Pi's optional.  This is inefficient, though, and precludes
% use of the matlab compiler.  Without them, we use feval and the number of Pi's must be
% changed with the editor for each application.  Places where this is required are marked
% with ARGLIST comments
%for i=nargin-6:-1:1
%   tailstr=[ ',P' num2str(i)  tailstr];
%   stailstr=[' P' num2str(i) stailstr];
%end
f0 = feval(fcn,x0,varargin{:});
%ARGLIST
%f0 = feval(fcn,x0,P1,P2,P3,P4,P5,P6,P7,P8,P9,P10,P11,P12,P13);
% disp('first fcn in csminwel.m ----------------') % Jinill on 9/5/95
if f0>1e50
  fprintf(fid+(fid==1),'Bad initial parameter.\n');
  return
else
  if verbose==1, fprintf(fid,'Iteration %5i: f = %20.8f\n',itct,f0); end
end
if NumGrad
  if length(grad)==0
    [g badg] = numgradvb(fcn,x0,verbose,fid,varargin{:});
    %ARGLIST
    %[g badg] = numgrad(fcn,x0,P1,P2,P3,P4,P5,P6,P7,P8,P9,P10,P11,P12,P13);
  else
    badg=any(find(grad==0));
    g=grad;
  end
  %numgrad(fcn,x0,P1,P2,P3,P4);
else
  [g badg] = feval(grad,x0,varargin{:});
  %ARGLIST
  %[g badg] = feval(grad,x0,P1,P2,P3,P4,P5,P6,P7,P8,P9,P10,P11,P12,P13);
end
retcode3=101;
x=x0;
f=f0;
H=H0;
cliff=0;
while ~done
  g1=[]; g2=[]; g3=[];
  %addition fj. 7/6/94 for control
  if verbose==2
    fprintf(fid,'-----------------\n');
    % vb = ShowText(vb,'f and x at the beginning of new iteration\n');
    fprintf(fid,'f at the beginning of new iteration, %20.10f\n',f);
    fprintf(fid,'x = %15.8g %15.8g %15.8g %15.8g\n',x);
    fprintf(fid,'\n');
  end
  %-------------------------
  itct=itct+1;
  [f1 x1 fc retcode1] = csminitvb(fcn,x,f,g,badg,H,verbose,fid,varargin{:});
  %ARGLIST
  %[f1 x1 fc retcode1] = csminit(fcn,x,f,g,badg,H,P1,P2,P3,P4,P5,P6,P7,...
  %           P8,P9,P10,P11,P12,P13);
  % itct=itct+1;
  fcount = fcount+fc;
  % erased on 8/4/94
  % if (retcode == 1) | (abs(f1-f) < crit)
  %    done=1;
  % end
  % if itct > nit
  %    done = 1;
  %    retcode = -retcode;
  % end
  if retcode1 ~= 1
    if retcode1==2 || retcode1==4
      wall1=1; badg1=1;
    else
      if NumGrad
        [g1 badg1] = numgradvb(fcn,x1,verbose,fid,varargin{:});
        %ARGLIST
        %[g1 badg1] = numgrad(fcn, x1,P1,P2,P3,P4,P5,P6,P7,P8,P9,...
        %                P10,P11,P12,P13);
      else
        [g1 badg1] = feval(grad,x1,varargin{:});
        %ARGLIST
        %[g1 badg1] = feval(grad, x1,P1,P2,P3,P4,P5,P6,P7,P8,P9,...
        %                P10,P11,P12,P13);
      end
      wall1=badg1;
      % g1
      %          save g1 g1 x1 f1 varargin;
      %ARGLIST
      %save g1 g1 x1 f1 P1 P2 P3 P4 P5 P6 P7 P8 P9 P10 P11 P12 P13;
    end
    if wall1 && (length(H) > 1)%
      % Bad gradient or back and forth on step length.  Possibly at
      % cliff edge.  Try perturbing search direction if problem not 1D
      %
      %fcliff=fh;xcliff=xh;
      Hcliff=H+diag(diag(H).*rand(nx,1));
      if verbose>0,fprintf(fid+(fid==1),'Cliff.  Perturbing search direction.\n');end
      [f2 x2 fc retcode2] = csminitvb(fcn,x,f,g,badg,Hcliff,verbose,fid,varargin{:});
      %ARGLIST
      %[f2 x2 fc retcode2] = csminit(fcn,x,f,g,badg,Hcliff,P1,P2,P3,P4,...
      %     P5,P6,P7,P8,P9,P10,P11,P12,P13);
      fcount = fcount+fc; % put by Jinill
      if  f2 < f
        if retcode2==2 || retcode2==4
          wall2=1; badg2=1;
        else
          if NumGrad
            [g2 badg2] = numgradvb(fcn,x2,verbose,fid,varargin{:});
            %ARGLIST
            %[g2 badg2] = numgrad(fcn, x2,P1,P2,P3,P4,P5,P6,P7,P8,...
            %      P9,P10,P11,P12,P13);
          else
            [g2 badg2] = feval(grad,x2,varargin{:});
            %ARGLIST
            %[g2 badg2] = feval(grad,x2,P1,P2,P3,P4,P5,P6,P7,P8,...
            %      P9,P10,P11,P12,P13);
          end
          wall2=badg2;
          % g2
          %                badg2
          %                save g2 g2 x2 f2 varargin
          %ARGLIST
          %save g2 g2 x2 f2 P1 P2 P3 P4 P5 P6 P7 P8 P9 P10 P11 P12 P13;
        end
        if wall2
          if verbose>0,fprintf(fid+(fid==1),'Cliff again. Try traversing\n');end
          if norm(x2-x1) < 1e-13
            f3=f; x3=x; badg3=1;retcode3=101;
          else
            gcliff=((f2-f1)/((norm(x2-x1))^2))*(x2-x1);
            if(size(x0,2)>1), gcliff=gcliff'; end
            [f3 x3 fc retcode3] = csminitvb(fcn,x,f,gcliff,0,eye(nx),verbose,fid,varargin{:});
            %ARGLIST
            %[f3 x3 fc retcode3] = csminit(fcn,x,f,gcliff,0,eye(nx),P1,P2,P3,...
            %         P4,P5,P6,P7,P8,...
            %      P9,P10,P11,P12,P13);
            fcount = fcount+fc; % put by Jinill
            if retcode3==2 || retcode3==4
              wall3=1; badg3=1;
            else
              if NumGrad
                [g3 badg3] = numgradvb(fcn, x3,verbose,fid,varargin{:});
                %ARGLIST
                %[g3 badg3] = numgrad(fcn, x3,P1,P2,P3,P4,P5,P6,P7,P8,...
                %                        P9,P10,P11,P12,P13);
              else
                [g3 badg3] = feval(grad,x3,varargin{:});
                %ARGLIST
                %[g3 badg3] = feval(grad,x3,P1,P2,P3,P4,P5,P6,P7,P8,...
                %                         P9,P10,P11,P12,P13);
              end
              wall3=badg3;
              % g3
              %                      badg3
              %                      save g3 g3 x3 f3 varargin;
              %ARGLIST
              %save g3 g3 x3 f3 P1 P2 P3 P4 P5 P6 P7 P8 P9 P10 P11 P12 P13;
            end
          end
        else
          f3=f; x3=x; badg3=1; retcode3=101;
        end
      else
        f3=f; x3=x; badg3=1;retcode3=101;
      end
    else
      % normal iteration, no walls, or else 1D, or else we're finished here.
      f2=f; f3=f; badg2=1; badg3=1; retcode2=101; retcode3=101;
    end
  else
    f2=f;f3=f;f1=f;retcode2=retcode1;retcode3=retcode1;
  end
  %how to pick gh and xh
  if f3 < f - crit && badg3==0
    ih=3;
    fh=f3;xh=x3;gh=g3;badgh=badg3;retcodeh=retcode3;
  elseif f2 < f - crit && badg2==0
    ih=2;
    fh=f2;xh=x2;gh=g2;badgh=badg2;retcodeh=retcode2;
  elseif f1 < f - crit && badg1==0
    ih=1;
    fh=f1;xh=x1;gh=g1;badgh=badg1;retcodeh=retcode1;
  else
    [fh,ih] = min([f1,f2,f3]);
    if verbose==2,fprintf(fid,'ih = %d',ih);end
    %eval(['xh=x' num2str(ih) ';'])
    switch ih
      case 1
        xh=x1;
      case 2
        xh=x2;
      case 3
        xh=x3;
    end %case
    %eval(['gh=g' num2str(ih) ';'])
    %eval(['retcodeh=retcode' num2str(ih) ';'])
    retcodei=[retcode1,retcode2,retcode3];
    retcodeh=retcodei(ih);
    if exist('gh')
      nogh=isempty(gh);
    else
      nogh=1;
    end
    if nogh
      if NumGrad
        [gh badgh] = numgradvb(fcn,xh,verbose,fid,varargin{:});
      else
        [gh badgh] = feval(grad, xh,varargin{:});
      end
    end
    badgh=1;
  end
  %end of picking
  %ih
  %fh
  %xh
  %gh
  %badgh
  stuck = (abs(fh-f) < crit);
  if (~badg)&&(~badgh)&&(~stuck)
    H = bfgsi(H,gh-g,xh-x);
  end
  if verbose==1
    fprintf(fid,'Iteration %5i: f = %20.8f, change = %20.8f\n',itct,fh,fh-f);
  elseif verbose==2
    fprintf(fid,'----\n');
    fprintf(fid,'Improvement on iteration %d = %20.8f\n\n',itct,f-fh);
  end
  % if Verbose
  if itct > nit
    if verbose>0,fprintf(fid+(fid==1),'iteration count termination\n');end
    done = 1;
  elseif stuck
    if verbose>0,fprintf(fid+(fid==1),'improvement < crit termination\n');end
    done = 1;
  end
  rc=retcodeh;
  if rc == 1
    if verbose>0,fprintf(fid+(fid==1),'zero gradient\n');end
  elseif rc == 6
    if verbose>0,fprintf(fid+(fid==1),'smallest step still improving too slow, reversed gradient\n');end
  elseif rc == 5
    if verbose>0,fprintf(fid+(fid==1),'largest step still improving too fast\n');end
  elseif (rc == 4) || (rc==2)
    if verbose>0,fprintf(fid+(fid==1),'back and forth on step length never finished\n');end
  elseif rc == 3
    if verbose>0,fprintf(fid+(fid==1),'smallest step still improving too slow\n');end
  elseif rc == 7
    if verbose>0,fprintf(fid+(fid==1),'warning: possible inaccuracy in H matrix\n');end
  end
  % end
  f=fh;
  x=xh;
  g=gh;
  badg=badgh;
end
% what about making an m-file of 10 lines including numgrad.m
% since it appears three times in csminwel.m

if ischar(fn),fclose(fid);end

end % function

%% -----------------------------------------------------------------------------

%% Function: csminitvb

function [fhat,xhat,fcount,retcode] = csminitvb(fcn,x0,f0,g0,badg,H0,verbose,fid,varargin)
% [fhat,xhat,fcount,retcode] = csminit(fcn,x0,f0,g0,badg,H0,...
%                                       P1,P2,P3,P4,P5,P6,P7,P8)
% retcodes: 0, normal step.  5, largest step still improves too fast.
% 4,2 back and forth adjustment of stepsize didn't finish.  3, smallest
% stepsize still improves too slow.  6, no improvement found.  1, zero
% gradient.
%---------------------
% Modified 7/22/96 to omit variable-length P list, for efficiency and compilation.
% Places where the number of P's need to be altered or the code could be returned to
% its old form are marked with ARGLIST comments.
%
% Fixed 7/17/93 to use inverse-hessian instead of hessian itself in bfgs
% update.
%
% Fixed 7/19/93 to flip eigenvalues of H to get better performance when
% it's not psd.
%
%tailstr = ')';
%for i=nargin-6:-1:1
%   tailstr=[ ',P' num2str(i)  tailstr];
%end
%ANGLE = .03;
ANGLE = .005;
%THETA = .03;
THETA = .3; %(0<THETA<.5) THETA near .5 makes long line searches, possibly fewer iterations.
FCHANGE = 1000;
MINLAMB = 1e-9;
% fixed 7/15/94
% MINDX = .0001;
% MINDX = 1e-6;
MINDFAC = .01;
fcount=0;
lambda=1;
xhat=x0;
f=f0;
fhat=f0;
g = g0;
gnorm = norm(g);
%
if (gnorm < 1.e-12) && ~badg % put ~badg 8/4/94
   retcode =1;
   dxnorm=0;
   % gradient convergence
else
   % with badg true, we don't try to match rate of improvement to directional
   % derivative.  We're satisfied just to get some improvement in f.
   %
   %if(badg)
   %   dx = -g*FCHANGE/(gnorm*gnorm);
   %  dxnorm = norm(dx);
   %  if dxnorm > 1e12
   %     disp('Bad, small gradient problem.')
   %     dx = dx*FCHANGE/dxnorm;
   %   end
   %else
   % Gauss-Newton step;
   %---------- Start of 7/19/93 mod ---------------
   %[v d] = eig(H0);
   %toc
   %d=max(1e-10,abs(diag(d)));
   %d=abs(diag(d));
   %dx = -(v.*(ones(size(v,1),1)*d'))*(v'*g);
   dx = -H0*g;
   dxnorm = norm(dx);
   if dxnorm > 1e12
      if verbose>0, fprintf(fid+(fid==1),'Near-singular H problem.\n');end
      dx = dx*FCHANGE/dxnorm;
   end
   dfhat = dx'*g0;
   %end
   if ~badg
      % test for alignment of dx with gradient and fix if necessary
      a = -dfhat/(gnorm*dxnorm);
      if a<ANGLE
         dx = dx - (ANGLE*dxnorm/gnorm+dfhat/(gnorm*gnorm))*g;
         dfhat = dx'*g;
         dxnorm = norm(dx);
         if verbose>0,fprintf(fid+(fid==1),'Correct for low angle: %g\n',a);end
      end
   end
   if verbose==2,fprintf(fid,'Predicted improvement: %20.8f\n',-dfhat/2);end
   %
   % Have OK dx, now adjust length of step (lambda) until min and
   % max improvement rate criteria are met.
   done=0;
   factor=3;
   shrink=1;
   lambdaMin=0;
   lambdaMax=inf;
   lambdaPeak=0;
   fPeak=f0;
   lambdahat=0;
   while ~done
      if size(x0,2)>1
         dxtest=x0+dx'*lambda;
      else
         dxtest=x0+dx*lambda;
      end
      % home
      f = feval(fcn,dxtest,varargin{:});
      %ARGLIST
      %f = feval(fcn,dxtest,P1,P2,P3,P4,P5,P6,P7,P8,P9,P10,P11,P12,P13);
      % f = feval(fcn,x0+dx*lambda,P1,P2,P3,P4,P5,P6,P7,P8);
      % conditional on Verbose added by Vasco Curdia
      if verbose==2,fprintf(fid,'lambda = %10.5g; f = %20.8f\n',lambda,f);end
      %debug
      %disp(sprintf('Improvement too great? f0-f: %g, criterion: %g',f0-f,-(1-THETA)*dfhat*lambda))
      if f<fhat
         fhat=f;
         xhat=dxtest;
         lambdahat = lambda;
      end
      fcount=fcount+1;
      shrinkSignal = (~badg & (f0-f < max([-THETA*dfhat*lambda 0]))) | (badg & (f0-f) < 0) ;
      growSignal = ~badg & ( (lambda > 0)  &  (f0-f > -(1-THETA)*dfhat*lambda) );
      if  shrinkSignal  &&   ( (lambda>lambdaPeak) || (lambda<0) )
         if (lambda>0) && ((~shrink) || (lambda/factor <= lambdaPeak))
            shrink=1;
            factor=factor^.6;
            while lambda/factor <= lambdaPeak
               factor=factor^.6;
            end
            %if (abs(lambda)*(factor-1)*dxnorm < MINDX) | (abs(lambda)*(factor-1) < MINLAMB)
            if abs(factor-1)<MINDFAC
               if abs(lambda)<4
                  retcode=2;
               else
                  retcode=7;
               end
               done=1;
            end
         end
         if (lambda<lambdaMax) && (lambda>lambdaPeak)
            lambdaMax=lambda;
         end
         lambda=lambda/factor;
         if abs(lambda) < MINLAMB
            if (lambda > 0) && (f0 <= fhat)
               % try going against gradient, which may be inaccurate
               lambda = -lambda*factor^6;
            else
               if lambda < 0
                  retcode = 6;
               else
                  retcode = 3;
               end
               done = 1;
            end
         end
      elseif  (growSignal && lambda>0) ||  (shrinkSignal && ((lambda <= lambdaPeak) && (lambda>0)))
         if shrink
            shrink=0;
            factor = factor^.6;
            %if ( abs(lambda)*(factor-1)*dxnorm< MINDX ) | ( abs(lambda)*(factor-1)< MINLAMB)
            if abs(factor-1)<MINDFAC
               if abs(lambda)<4
                  retcode=4;
               else
                  retcode=7;
               end
               done=1;
            end
         end
         if ( f<fPeak ) && (lambda>0)
            fPeak=f;
            lambdaPeak=lambda;
            if lambdaMax<=lambdaPeak
               lambdaMax=lambdaPeak*factor*factor;
            end
         end
         lambda=lambda*factor;
         if abs(lambda) > 1e20;
            retcode = 5;
            done =1;
         end
      else
         done=1;
         if factor < 1.2
            retcode=7;
         else
            retcode=0;
         end
      end
   end
end
if verbose==2,fprintf(fid,'Norm of dx %10.5g\n', dxnorm);end

end % function csminitvb


%% -----------------------------------------------------------------------------


