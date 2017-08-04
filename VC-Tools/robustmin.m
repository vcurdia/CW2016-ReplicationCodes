function Out = robustmin(fcn,x0,varargin)

% robustmin
%
% This routine builds on Chris Sims routine csminwel to make its results
% more robust. The actual implementation uses csminwelvb.
%
% More in detail the robustness comes from running csminwel several times
% in succession, in which each guess will use the previous output, shock it
% and use the result as the new guess value. Furthermore the initial
% inverse hessian for each iteration is the inverse hessian from the
% previous one.
%
% Usage:
%   Out = robustmin(fcn,x0)
%   Out = robustmin(fcn,x0,op)
%
% Inputs:
%   fcn
%     name or handle of function to be minimized
%   x0
%     Column vector with the guess values
%   varargin (OPTIONAL):
%     Can be set in one of two ways:
%     1) Cell array with input arguments for the function being minimized.
%        No other minimization parameters are set.
%     2) Structure with new values for the parameters of the minimization
%        and robustness check. The structure fields need to be among the
%        following options, but do not have to include all of them - this
%        way in a specific application only needs to set the parameters for
%        which the default is not good.
%        2.1) parameters to be used in csminwelvb
%          op.H0
%            Initial hessian. Default: eye(length(x0))
%          op.grad
%            Name of gradient function. Default: []
%          op.crit
%            Criterion for convergence. Default: 1e-7
%          op.nit
%            Maximum number of iterations. Default: 1000
%          op.verbose
%            Choose degree of verbosity. 0 does not show anything. 1 shows only
%            function value at each iteration and other warning messages. 2
%            shows everything. Default: 2
%          op.varargin [cell array]
%            Additional arguments required to evaluate function. This is
%            required if the function requires arguments. Default: {}
%        2.2) parameters needed for robustness part
%          op.Rcrit
%            Criterion for robustness convergence. Default: 1e-7
%          op.Ritmax
%            Maximum number of robustness iterations. If set to zero, no
%            robustness is performed. Default: 30
%          op.Ritmin
%            Minimum number of robustness iterations. Default: 10
%          op.Ritnoimp
%            Maximum number of iterations without improvement before
%            robustness is declared. Default: 5
%          op.RH0
%            If set to 1 then the robust shock is based on the initial H0,
%            if set to 0 then the robust shock is based on the latest H.
%            Default: 1
%          op.Rscaledown
%            Scaling of previous inverse hessian used for drawing candidate
%            guess value. Default: 0.25 (equivalent to 50% of SE)
%        2.3) Other parameters
%          op.AltHessian
%            If set to 1, Hessian is computed using vcHessian,
%            if set to 0, Hessian is computed in csminwel.
%            Default: 0
% Output:
%   The output of the minimization is introduced in a structure, as follows
%   1) Output typical of vccsminwel
%     Out.x
%       Argmin for this function
%     Out.f
%       Value of function at its minimum
%     Out.g
%       Gradient at minimum
%     Out.H
%       Hessian at minimium
%     Out.itct
%       Iteration count from last call of csminwel at minimum
%     Out.fc
%       Number of gradients computed in csminit at minimum
%     Out.rc
%       Return code from csminwel at minimum
%     Out.rcMsg
%       Message explaining return code from csminwel
%       rc = 0 => Normal solution
%       rc = 1 => Zero gradient
%       rc = 2 or rc = 4 => Back and forth on step length never finished
%       rc = 3 => Smallest step still improving too slow
%       rc = 5 => Largest step still improving too fast
%       rc = 6 => Smallest step still improving too slow, reversed gradient
%       rc = 7 => Warning: possible inaccuracy in H matrix
%       rc = 8 => Invalid initial guess
%     Out.fn
%       If fn specified it is the filename with verbose.
%   2) Output from robust analysis
%     Out.Ritct
%       Number of iterations of robustness
%     Out.Rrc
%       Return code for robustness
%     Out.RrcMsg
%       Message describing return code for robustness
%       Rrc = 0 => no robustness
%       Rrc = 1 => improvement in function below Rcrit for Ritnoimp
%                  consecutive iterations
%       Rrc = 2 => maximum number of iterations exceeded
%   3) Include the csminwel results for initial and robust iterations
%     Out.MinOutput
%       Structure array with 1+Ritct elements. The fields are the output of
%       csminwel for the corresponding minimization.
%   3) also list parameters used in minimization
%     Out.op
%       Structure with all the parameters used, including guess value
%
%
% See also:
% csminwel, vccsminwel csminit, numgrad, csolve, vcnumjacobian, vchessian
%
% ..............................................................................
%
% Created: February 5, 2008 by Vasco Curdia
%
% Copyright 2008-2017 by Vasco Curdia

%% -----------------------------------------------------------------------------

%% Default parameter values
% parameters needed for vccsminwel
op.x0 = x0;
op.H0 = eye(length(x0));
op.grad = [];
op.crit = 1e-7;
op.nit = 1000;
op.verbose = 1;
op.varargin = {};

% parameters needed for robustness part
op.Rcrit = 1e-7;
op.Ritmax = 30;
op.Ritmin = 10;
op.Ritnoimp = 5;
op.RH0 = 1;
op.Rscaledown = 0.25;

% other parameters
op.AltHessian = 0;
op.MatFn = [];
op.LogFn = 1;

op = updateoptions(op,varargin{:});

Out.op = op;

%% Check consistency
% op.Ritmin = min(op.Ritmin,op.Ritmax);
% op.Ritnoimp = min(op.Ritnoimp,op.Ritmax);

%% Check output destination
if ischar(op.LogFn)
    fid = fopen(op.LogFn,'wt');
else
    fid = op.LogFn;
end

%% Check if need to save results and where
isSaveMat = ~isempty(op.MatFn);

%% Save op in output structure
Out.op = op;

%% Display parameters used:
fprintf(fid,'\nGuess vector:\n');
fprintf(fid,'%20.10f\n',x0);
fprintf(fid,'\n');
Out.x0 = x0;

%% -----------------------------------------------------------------------------

%% Check if it is good guess
f0 = feval(fcn,x0,op.varargin{:});
Out.f0=f0;
if isSaveMat,save(op.MatFn,'Out');end
if f0>1e50
    fprintf(fid+(fid==1),'\nWARNING: Invalid initial guess. Cannot proceed.\n');
    Out.x = x0;
    Out.f = f0;
    Out.g = [];
    Out.H = [];
    Out.itctr = 0;
    Out.fcr = 0;
    Out.rc = 8;
    Out.rcMsg = 'Invalid initial guess';
    Out.Ritct = 0;
    Out.Rrc = 0;
    Out.RrcMsg = 'No Robustness performed';
    Out.MinOutput = [];
    Out.op = op;
    if isSaveMat,save(op.MatFn,'Out');end
    return
end

%% Make first minimization
MinOutputFields = {'fh','xh','gh','Hh','itct','fc','rc'};
nMinOutputFields = length(MinOutputFields);
[fh,xh,gh,Hh,itct,fc,rc] = ...
    csminwelvb(fcn,x0,op.H0,op.grad,op.crit,op.nit,op.verbose,fid,...
               op.varargin{:});
if op.AltHessian
    Hh = inv(vcHessian(fcn,xh,op.varargin{:}));
end
for j=1:nMinOutputFields
    MinOutput(1).(MinOutputFields{j}) = eval(MinOutputFields{j});
end
MinOutput(1).rcMsg = rcErrorMsg(rc);
xhr = xh;
fhr = fh; fh0 = fh;
Hhr = Hh;
ghr = gh;
Hh_backup = op.H0;
x0_backup = x0;
rcr = rc;
itctr = itct;
fcr = fc;
Out.x = xhr;
Out.f = fhr;
Out.H = Hhr;
Out.g = gh;
Out.itct = itctr;
Out.fc = fcr;
Out.rc = rcr;
Out.rcMsg = rcErrorMsg(rcr);
if exist('fn','var'),Out.fn = fn;end
Out.MinOutput = MinOutput;
if isSaveMat,save(op.MatFn,'Out');end

%% -----------------------------------------------------------------------------

%% Robustness
noimp = 0;
for rit=1:op.Ritmax
    fprintf(fid,'\nRobust Iteration %.0f\n',rit);
% check for semi positive definite inverse hessian
    [Ur,Dr] = eig((Hhr+Hhr')/2);
    Dr = diag(Dr);
    tol = eps(max(Dr)) * length(Dr);
    idxD = (abs(Dr) > tol);
    Dr = Dr(idxD);
    negeig = sum(Dr<0); % number of negative eigenvalues
    if negeig~=0
        fprintf(fid+(fid==1),'Warning: Previous Hessian is not positive semi definite!\n');
        fprintf(fid+(fid==1),'         Using earlier best alternative...\n');
        Hhr = Hh_backup;
    else
        Hh_backup = Hhr;
    end
% choose H for robust shock
    if op.RH0
        RShockVar = op.H0;
    else
        RShockVar = Hhr;
    end
    x0r = mvnrnd(xhr,op.Rscaledown*RShockVar)';
    
% check whether it is an acceptable candidate
    GoodC = (feval(fcn,x0r,op.varargin{:})<1e50);
    for j=1:1000
        if GoodC, break, end
        x0r = mvnrnd(xhr,op.Rscaledown*RShockVar)';
        GoodC = (feval(fcn,x0r,op.varargin{:})<1e50);
    end
    if ~GoodC
        fprintf(fid+(fid==1),'Warning: could not find suitable candidate after 1000 attempts.\n');
        fprintf(fid+(fid==1),'         Using last good candidate...\n');
        x0r = x0_backup;
    end
    
    [fh,xh,gh,Hh,itct,fc,rc] = ...
        csminwelvb(fcn,x0r,Hhr,op.grad,op.crit,op.nit,op.verbose,fid,...
                   op.varargin{:});
    if op.AltHessian
        Hh = inv(vcHessian(fcn,xh,op.varargin{:}));
    end
    for j=1:nMinOutputFields
        MinOutput(1+rit).(MinOutputFields{j}) = eval(MinOutputFields{j});
    end
    MinOutput(1+rit).rcMsg = rcErrorMsg(rc);
  
    fprintf(fid,'\nResults from robustness iteration %.0f\n',rit);
    fprintf(fid,'old function value: %.10f\n', fhr);
    fprintf(fid,'new function value: %.10f\n', fh);
    dfh=fh-fhr;
    fprintf(fid,'change in function value: %.10f\n\n', dfh);
    rf(rit) = fh;
    df(rit) = dfh;
    if -dfh>=op.Rcrit
        xhr = xh;
        fhr = fh;
        Hhr = Hh;
        ghr = gh;
        rcr = rc;
        itctr = itct;
        fcr = fc;
        noimp=0;
    else
        noimp = noimp + 1;
        if (noimp>=op.Ritnoimp) && (rit>=op.Ritmin)
            Rrc = 1;
            RrcMsg = sprintf('Improvement in function below %g for %.0f consecutive iterations',...
                             op.Rcrit, op.Ritnoimp);
            break
        end
    end
    Out.x = xhr;
    Out.f = fhr;
    Out.H = Hhr;
    Out.g = gh;
    Out.itct = itctr;
    Out.fc = fcr;
    Out.rc = rcr;
    Out.rcMsg = rcErrorMsg(rcr);
    if exist('fn','var'),Out.fn = fn;end
    Out.MinOutput = MinOutput;
    if isSaveMat,save(op.MatFn,'Out');end
end

%% Checks
if op.Ritmax==0
    Rrc = 0;
    RrcMsg = 'No robustness performed';
    rf = [];
    rit = 0;
elseif rit==op.Ritmax
    Rrc = 2;
    RrcMsg = 'Maximum number of iterations exceeded';
end

%% print history of posterior values and improvements
% % fprintf('\nFunction value before robustness: %f\n', fh0)
% fprintf('Initial iterarion, function value: %.10f',fh0)
% fprintf('(Stopped at iteration %.0f)\n',MinOutput(1).itct)
% if Rrc~=0
%     fprintf('\nResults from robustness check:')
%     fprintf('\n------------------------------\n\n')
%     for rit=1:length(rf)
%         fprintf('Robustness iterarion %.0f, function value: %.10f, change %.10f',rit,rf(rit),df(rit))
%         fprintf('(Stopped at iteration %.0f)\n',MinOutput(1+rit).itct)
%     end
%     fprintf('\nFunction value after robustness: %f\n', fhr)
% end
% fprintf('\n%s\n',RrcMsg)

fprintf(fid,'\nRobustness analysis for minimization');
fprintf(fid,'\n------------------------------------\n');
fprintf(fid,'Iteration %3.0f: function value: %15.8f',0,MinOutput(1).fh);
% fprintf('         %15s','')
fprintf(fid,' change: %15.8f',MinOutput(1).fh-f0);
fprintf(fid,' stopped at iteration %.0f\n',MinOutput(1).itct);
itbest.fh = MinOutput(1).fh;
itbest.idx = 1;
for jr=2:length(MinOutput)
    fprintf(fid,'Iteration %3.0f: function value: %15.8f',jr-1,MinOutput(jr).fh);
    itchg = MinOutput(jr).fh-itbest.fh;
    fprintf(fid,' change: %15.8f',itchg);
    fprintf(fid,' stopped at iteration %.0f\n',MinOutput(jr).itct);
    if itchg<0
        itbest.fh = MinOutput(jr).fh;
        itbest.idx = jr;
    end
end
fprintf(fid,'\nBest iteration: %.0f',itbest.idx-1);
fprintf(fid,'\nBest iteration function value: %.8f\n',itbest.fh);
fprintf(fid,'\nRobustness change over initial min: %.8f',itbest.fh-MinOutput(1).fh);
fprintf(fid,'\nRobustness message: %s\n',RrcMsg);
fprintf(fid,'\nBest iteration message: %s\n\n',MinOutput(itbest.idx).rcMsg);

%% check for positive semi definite hessian
[U,D] = eig((Hhr+Hhr')/2);
D = diag(D);
tol = eps(max(D)) * length(D);
idxD = (abs(D) > tol);
D = D(idxD);
negeig = sum(D<0); % number of negative eigenvalues
if negeig~=0
    fprintf(fid+(fid==1),'Warning: Hessian is not positive semi definite!\n');
    fprintf(fid+(fid==1),'         Using last good one...\n');
    Hhr = Hh_backup;
    if all(all(Hh_backup==op.H0))
        fprintf(fid+(fid==1),'Backup Hessian is initial guess!\n');
    end
end

%% -----------------------------------------------------------------------------

%% Prepare output

% Output typical of vccsminwel, for the optimal solution
Out.x = xhr;
Out.f = fhr;
Out.H = Hhr;
Out.g = gh;
Out.itct = itctr;
Out.fc = fcr;
Out.rc = rcr;
Out.rcMsg = rcErrorMsg(rcr);
if exist('fn','var'),Out.fn = fn;end

% Output from robust analysis
if Rrc~=0
    Out.Ritct = rit;
else
    Out.Ritct = 0;
end
Out.Rrc = Rrc;
Out.RrcMsg = RrcMsg;

% for each minimization save the output:
Out.MinOutput = MinOutput;

if ischar(op.LogFn), fclose(fid); end

%% -----------------------------------------------------------------------------

%% Subfunction: rcErrorMsg(rc)
function rcMsg = rcErrorMsg(rc)
if rc == 1
    rcMsg = 'Zero gradient';
elseif rc == 6
    rcMsg = 'Smallest step still improving too slow, reversed gradient';
elseif rc == 5
    rcMsg = 'Largest step still improving too fast';
elseif (rc == 4) || (rc==2)
    rcMsg = 'Back and forth on step length never finished';
elseif rc == 3
    rcMsg = 'Smallest step still improving too slow';
elseif rc == 7
    rcMsg = 'Warning: possible inaccuracy in H matrix';
else
    rcMsg = 'Normal solution';
end

%% -----------------------------------------------------------------------------



