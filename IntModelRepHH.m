function IntModelRepHH(varargin)

% IntModelRepHH
%
% Solves for LQ optimal policy in Curdia and Woodford (2008).
% This version is the RepHH model, nested in the full model.
%
% Usage:
%   IntModelRepHH
%   IntModelRepHH(...,OptionName,...)
%
% Convention for refering to variables in equations:
%   x_t   refers to x{t}
%   x_tF  refers to x{t+1}
%   x_tL refers to x{t-1}
%   x_ss  refers to the steady state level of 'x'
%   (where x refers to some variable with name 'x')
%
% Required m-files:
%   - symbolic toolbox
%   - LQGenSymVar
%   - LQ
%   - LQSolveREE
%   - LQCheckSOC
%   - gensys.m, available in Chris Sims's website
%   - csolve.m, available in Chris Sims's website
%
% See also:
% LQ, LQSolveREE, LQCheckSOC, LQCheckSOCold, LQGenSymVar  
%
% .........................................................................
%
% Created: January 24, 2008
% Updated: July 29, 2010
% by Vasco Curdia and Michael Woodford

%% ------------------------------------------------------------------------

%% preamble
% tic
ttic = toc();

nsteps = 25;
NumPrecision = 1e-10;
SolveIterations = 1000;

%% Defaults
isPers = 1;
isNatVars = 0;
isEnd = 1;
isNoRes = 0;
isNoSpread = 0;
isNoDist = 0;
isGDebt = 0;
isSmSigma = 0;
isLowSigma = 0;

%% check arguments
if ismember('NoPers',varargin), isPers = 0; end
if ismember('Pers',varargin), isPers = 1; end
if ismember('NatVars',varargin), isNatVars = 1; end
if ismember('NoNatVars',varargin), isNatVars = 0; end
if ismember('Exo',varargin), isEnd = 0; end
if ismember('End',varargin), isEnd = 1; end
if ismember('NoRes',varargin), isNoRes = 1; end
if ismember('NoSpread',varargin), isNoSpread = 1; end
if ismember('NoDist',varargin), isNoDist = 1; end
if ismember('GDebt',varargin), isGDebt = 1; end
if ismember('SmSigma',varargin), isSmSigma = 1; end
if ismember('LowSigma',varargin), isLowSigma = 1; end

%% Exercise type
ExerciseName = 'RepHH';
if isPers
    ExerciseName = [ExerciseName,'_Pers'];
else
    ExerciseName = [ExerciseName,'_NoPers'];
end
if isEnd
    ExerciseName = [ExerciseName,'_End'];
else
    ExerciseName = [ExerciseName,'_Exo'];
end
if isNoRes
    ExerciseName = [ExerciseName,'_NoRes'];
end
if isNoSpread
    ExerciseName = [ExerciseName,'_NoSpread'];
end
if isNoDist
    ExerciseName = [ExerciseName,'_NoDist'];
end
if isGDebt
    ExerciseName = [ExerciseName,'_GDebt'];
end
if isSmSigma
    ExerciseName = [ExerciseName,'_SmSigma'];
end
if isLowSigma
    ExerciseName = [ExerciseName,'_LowSigma'];
end
% DiaryFileName = [ExerciseName,'.log'];
% if exist(DiaryFileName,'file')
%     movefile([ExerciseName,'.log'],[ExerciseName,'_Old.log'])
% end
% diary(DiaryFileName)
fprintf('\n****%s****',repmat('*',1,length(ExerciseName)))
fprintf('\n*   %s   *',ExerciseName)
fprintf('\n****%s****\n\n',repmat('*',1,length(ExerciseName)))

%% ------------------------------------------------------------------------

%% define parameters
rd_ss = 0.01;
sigmabari = isLowSigma*0.9+~isLowSigma*0.16; %0.16
sigmabar = 1/sigmabari;
phi = 1/0.75;
alpha = 0.66;
omega_y = 0.473;
nu = (omega_y+1)/phi-1;
mu_p = 1.15;
theta = 1+1/(mu_p-1);
Yss = 1;
s_c = 0.7;
mu_w_ss = 1;
tau_ss = ~isNoDist*0.2+isNoDist*(1-mu_p*mu_w_ss);
rho_bg = isGDebt*0.1;
Hbar_ss = 1;
Z_ss = 1;
rho_csi = isPers*0.9;
rho_m = isPers*0.6;
phi_pi = 2; %2
phi_y = 1; %1

%% Endogenous variables
y = {...
    'Rd', 1.010000;
    'Pi', 1;
    'Y', 0.414274;
    'lambda', 1.263350;
    'Delta', 1;
    'K', 1.208245;
    'F', 1.208245;
    };
if isNatVars
    y = cat(1,y,{...
        % natural variables
        'Rrdn', 1.010000;
        'Yn', 0.414274;
        'lambdan', 1.263350;
        });
end

%% Exogenous variables
csi = {'xi_i','hZ','hmu_w','htau','hG','hb_g','hHbar','hCbar'};
S = diag([rho_m,rho_csi*ones(1,length(csi)-1)]); 

%% Lagrange multipliers for F constraints
nF = 3+isNatVars*2; for jF=1:nF,FLM(jF,:)={sprintf('FLM%.0f',jF),0};end

%% Lagrange multipliers for G constraints
nG = 3+isNatVars*1; for jG=1:nG,GLM(jG,:)={sprintf('GLM%.0f',jG),0};end

%% Generate symbolic variables
LQGenSymVar

%% ------------------------------------------------------------------------

%% Auxiliary definitions, used in later sections

% parameter values and steady state values
Rd_ss = 1+rd_ss;
Pi_ss = 1;
Delta_ss = 1;
Y_ss = Yss;
beta = 1/Rd_ss;
lambda_ss = mu_p*phi*mu_w_ss*Hbar_ss^(-nu)*(Y_ss/Z_ss)^(1+omega_y)/((1-tau_ss)*Y_ss);
s_g = 1-s_c;
sigma = sigmabar/s_c;
Cbar_ss = s_c*Y_ss*lambda_ss^sigma;
ctil_ss = s_c*Y_ss;
K_ss = mu_p*phi*mu_w_ss*Hbar_ss^(-nu)*(Y_ss/Z_ss)^(1+omega_y)/(1-alpha*beta);
F_ss = lambda_ss*(1-tau_ss)*Y_ss/(1-alpha*beta);
Rrdn_ss = Rd_ss/Pi_ss;
Yn_ss = Y_ss;
lambdan_ss = lambda_ss;

% shocks
Z_t = Z_ss*exp(hZ_t);
mu_w_t = mu_w_ss*exp(hmu_w_t);
tau_t = tau_ss+tau_ss*htau_t;
G_t = s_g*Y_ss+Y_ss*hG_t;
b_g_t = rho_bg*Y_ss+Y_ss*hb_g_t;
Hbar_t = Hbar_ss*exp(hHbar_t);
Cbar_t = Cbar_ss*exp(hCbar_t);

% auxiliary definitions
ctil_t = Cbar_t*lambda_t^(-sigma);
if isNatVars
    ctiln_t = Cbar_t*lambdan_t^(-sigma);
end

%% Utility function
U = ctil_t^(1-1/sigma)*Cbar_t^(1/sigma)/(1-1/sigma)...
    -1/(1+nu)*Hbar_t^(-nu)*(Y_t/Z_t)^(1+omega_y)*Delta_t;

%% F constraints
F = vpa([...
    Cbar_t*lambda_t^(-sigma)+G_t-Y_t;
	alpha*Delta_tL*Pi_t^(theta*(1+omega_y))+...
        (1-alpha)*((1-alpha*Pi_t^(theta-1))/(1-alpha))^(theta*(1+omega_y)/(theta-1))-...
        Delta_t;
	(F_t/K_t)^((theta-1)/(1+omega_y*theta))-(1-alpha*Pi_t^(theta-1))/(1-alpha);
    ]);
if isNatVars
    F = vpa([F;
        % natural variable equations
        Cbar_t*lambdan_t^(-sigma)+G_t-Yn_t;
        mu_p*phi*mu_w_ss*Hbar_t^(-nu)*(Yn_t/Z_t)^(1+omega_y)-lambdan_t*(1-tau_ss)*Yn_t
        ]);
end

%% G constraints
G = vpa([...
	Rd_t*beta*lambda_tF/Pi_tF-lambda_t;
	mu_p*phi*mu_w_t*Hbar_t^(-nu)*(Y_t/Z_t)^(1+omega_y)+...
        alpha*beta*Pi_tF^(theta*(1+omega_y))*K_tF-K_t;
	lambda_t*(1-tau_t)*Y_t+alpha*beta*Pi_tF^(theta-1)*F_tF-F_t;
    ]);
if isNatVars
    G = vpa([G;
        % natural variable equations
        Rrdn_t*beta*lambdan_tF-lambdan_t;
        ]);
end

%% ------------------------------------------------------------------------

%% Compute steady state
fprintf('\nSolving for steady state...')

%% generate the steady state system
ssSys = vpa(jacobian(U,y_t).'+...
        (FLM_ss*jacobian(F,y_t)).'+beta*(FLM_ss*jacobian(F,y_tL)).'+...
        (GLM_ss*jacobian(G,y_t)).'+beta^(-1)*(GLM_ss*jacobian(G,y_tF)).');
ssSys = [ssSys; F];
ssSys = [ssSys; G];
% plug in the steady state values in symbolic form
ssSys = subs(ssSys,[y_t,y_tF,y_tL,csi_t],[y_ss,y_ss,y_ss,csi_ss]);

%% prepare variables
idxPi = find(ismember({y(:).name},'Pi'));
ssSys1 = eval(ssSys([1:idxPi-1,idxPi+1:ny]));
nSys = length(ssSys1);
x = FLM;
x(end+(1:nG)) = GLM;
nx = length(x);

%% generate function for csolve
SolveFileName = sprintf('ssSys%s',ExerciseName);
fid=fopen([SolveFileName,'.m'],'w');
fprintf(fid,'function f=%s(x) \n',SolveFileName);
fprintf(fid,'f = ones(size(x));\n');
fprintf(fid,'for j=1:size(x,2)  \n');
for j=1:nSys
    fprintf(fid,'%s_ss = x(%.0f,j);  \n',x(j).name,j);
%     if ismember(x(j).name,{y(:).name})
%         fprintf(fid,'if %s_ss<0, f(:,j) = inf; continue, end\n',x(j).name);
%     end
end
for j=1:nSys
    fprintf(fid,['f(',int2str(j),',j) = ',char(ssSys1(j)),';  \n']);
end
fprintf(fid,'end  \n');
fclose(fid);

%% Solve
[x1,rc] = csolve(SolveFileName,[x.guess]',[],NumPrecision,SolveIterations);
if rc~=0, error(['Solution of steady state system is not normal, rc = ', int2str(rc)]), end
% [(1:nx)' feval(SolverFileName,x1)] % check system solution
delete([SolveFileName,'.m'])

%% evaluate variables (need to keep different notation, so as not to erase
% symbolic expressions at this stage)
% y_ss_1 = x1(1:ny)';
% FLM_ss_1 = x1(ny+1:ny+nF)';
% GLM_ss_1 = x1(ny+nF+1:ny+nF+nG)';
for j=1:nx
    eval(sprintf('%s_ss = %.16f;',x(j).name,x1(j)))
end
Phi = 1-(theta-1)/theta*(1-tau_ss)/mu_w_ss;
xi = (1-alpha)/alpha*(1-alpha*beta)/(1+omega_y*theta);
kappa = xi*(omega_y+sigmabari);

ssSys = eval(ssSys);
if ~all(abs(ssSys)<1e-6)
    fprintf('\nWARNING: system solution is not precise\n')
    [(1:length(ssSys))' ssSys]
end

%% Present steady state
fprintf('\nSteady state results:')
fprintf('\n=====================\n\n')
for j=1:ny
    fprintf('%15s_ss = %12.6f\n',y(j).name,eval([y(j).name,'_ss']))
end
for j=1:nF
    fprintf('%15s_ss = %12.6f\n',FLM(j).name,eval([FLM(j).name,'_ss']))
end
for j=1:nG
    fprintf('%15s_ss = %12.6f\n',GLM(j).name,eval([GLM(j).name,'_ss']))
end
disp(' ')

%% present some ratios
fprintf('\nSome ratios:')
fprintf('\n============\n\n')
xdisp({'beta','s_g','rho_bg','s_c',...
    'sigmabar','sigma','Cbar_ss','ctil_ss',...
    'phi_pi','phi_y'...
    'Phi','xi','kappa'})
disp(' ')

%% ------------------------------------------------------------------------

%% LQ Optimal Policy
fprintf('\nSolving for LQ Optimal Policy...\n')

%% Run LQ function to get symbolic matrices
[LQmat,LQk_t]=LQ(y_t,csi_t,U,F,G,beta,S,y_ss,FLM_ss,GLM_ss);

%% Fill LQ matrices with steady state and parameter vectors
for j=1:nx
    eval([x(j).name,'_ss = x1(',int2str(j),');'])
end
y_ss = eval(y_ss);
FLM_ss = eval(FLM_ss);
GLM_ss = eval(GLM_ss);
matnames = fieldnames(LQmat);
for j=1:length(matnames)
    LQmat.(matnames{j}) = reshape(eval(LQmat.(matnames{j})(:)),size(LQmat.(matnames{j})));
end

%% Solve the REE
[REE.LQ,LQz_t] = LQSolveREE(LQmat,LQk_t);

%% Check for the SOC
SOC = LQCheckSOC(beta,S,LQmat,1e-6,20);

%% ------------------------------------------------------------------------

%% Alternative rules
fprintf('\nSolving for REE of alternative rules...\n')

%% generate symbolic variables for rules
for j=1:ny
    eval(sprintf('syms %1$s %1$sF %1$sL',char(LQz_t(j))))
end

%% aux definitions
omega_tau = tau_ss/(1-tau_ss);
q_pi = theta/kappa*(omega_y+sigmabari+Phi*(1-sigmabari));
q_y = omega_y+sigmabari+Phi*(1-sigmabari)-Phi/sigmabar*(1/s_c-1)/(omega_y+sigmabari);
omega_1 = 1/q_y*(omega_y+sigmabari+Phi*(1-sigmabari));
omega_2 = Phi/s_c*sigmabari...
    /((omega_y+sigmabari)^2+Phi*((1-sigmabari)*(omega_y+sigmabari)-(1/s_c-1)*sigmabari));
omega_3 = (1-Phi)...
    /(omega_y+sigmabari+Phi*(1-sigmabari-(1/s_c-1)*sigmabari/(omega_y+sigmabari)));
omega_4 = omega_tau...
    /(omega_y+sigmabari+Phi*(1-sigmabari-(1/s_c-1)*sigmabari/(omega_y+sigmabari)));
q_t = ((1+omega_y)*hZ_t+nu*hHbar_t)/omega_y;
q_tL = ((1+omega_y)*hZ_tL+nu*hHbar_tL)/omega_y;
g_t = s_c*hCbar_t+hG_t;
g_tL = s_c*hCbar_tL+hG_tL;
hYnBW_t = (sigmabari*g_t+omega_y*q_t-hmu_w_t-omega_tau*htau_t)/(omega_y+sigmabari);
hYnBW_tL = (sigmabari*g_tL+omega_y*q_tL-hmu_w_tL-omega_tau*htau_tL)/(omega_y+sigmabari);
hYsBW_t = omega_1*hYnBW_t-omega_2*hG_t+omega_3*hmu_w_t+omega_4*htau_t;
hYsBW_tL = omega_1*hYnBW_tL-omega_2*hG_tL+omega_3*hmu_w_tL+omega_4*htau_tL;
hyBW_t = hY_t-hYsBW_t;
hyBW_tL = hY_tL-hYsBW_tL;

%% Alternative rule: Taylor
Rule.Taylor = phi_pi*hPi_t+phi_y/4*hY_t+xi_i_t-hRd_t;

%% Alternative rule: PiStab
Rule.PiStab = hPi_t;

%% Alternative rule: FlexTarget
lambda_x = q_y/q_pi/kappa;
% Rule.FlexTarget = hPi_t+lambda_x*(hyBW_t-hyBW_tL);
Rule.FlexTarget = hPi_t+lambda_x*((hY_t-hYsBW_t)-(hY_tL-hYsBW_tL));

%% Alternative rule: TaylorRn
% Rule.TaylorRn = hRrdn_t+phi_pi*hPi_t+phi_y/4*hY_t+xi_i_t-hRd_t;

%% Alternative rule: TaylorYn
Rule.TaylorYn = phi_pi*hPi_t+phi_y/4*(hY_t-hYn_t)+xi_i_t-hRd_t;

%% Alternative rule: TaylorYnBW
% Rule.TaylorYnBW = phi_pi*hPi_t+phi_y/4*hyBW_t+xi_i_t-hRd_t;

%% Alternative rule: TaylorRnYn
% Rule.TaylorRnYn = hRrdn_t+phi_pi*hPi_t+phi_y/4*(hY_t-hYn_t)+xi_i_t-hRd_t;

%% Run rules
AltPol = fieldnames(Rule);
nAltPol = length(AltPol);
for jPol=1:nAltPol
    fprintf('  %s\n',AltPol{jPol})
    REE.(AltPol{jPol}) = LQAltRuleLinear(Rule.(AltPol{jPol}),LQmat,S,LQz_t);
end

%% ------------------------------------------------------------------------

%% Generate IRFs
fprintf('\nGenerating IRFs...\n')
Policy = fieldnames(REE);
nPol = length(Policy);
zz = {y(:).name,csi{:},'Rb','Rrd','Rrb','cb','cs','w','wb','ws'};
if isNatVars
    zz = cat(2,zz,{'Rrbn','cbn','csn','wn','wbn','wsn'});
end
for j=1:length(zz)
    eval(sprintf('[tf,idx%1$s] = ismember(''%1$s'',zz);',zz{j}))
end
ShockSize = ones(ncsi,1);
ShockSize(ismember(csi,'htau')) = 1/tau_ss;
ShockSize(ismember(csi,'xi_i')) = 1/4;
ShockSize(ismember(csi,'hb_g')) = 4;
ShockSize = diag(ShockSize);
for jPol=1:nPol
    irf = REE.(Policy{jPol}).Phi2*ShockSize;
    for j=2:nsteps
        irf(:,:,j) = REE.(Policy{jPol}).Phi1*irf(:,:,j-1);
    end
    % select only model variables, excluding multipliers and LQ artificial variables
    irf = irf(1:ny+ncsi,:,:);
    % add lending rate
    irf(idxRb,:,:) = irf(idxRd,:,:);
    % add real deposit rate
    irf(idxRrd,:,:) = irf(idxRd,:,:)-cat(3,irf(idxPi,:,2:end),NaN(1,ncsi));
    % add real lending rate
    irf(idxRrb,:,:) = irf(idxRb,:,:)-cat(3,irf(idxPi,:,2:end),NaN(1,ncsi));
    % add cb
    irf(idxcb,:,:) = irf(idxhCbar,:,:)-sigma*irf(idxlambda,:,:);
    % add cs
    irf(idxcs,:,:) = irf(idxhCbar,:,:)-sigma*irf(idxlambda,:,:);
    % add w
    irf(idxw,:,:) = irf(idxhmu_w,:,:)-nu*irf(idxhHbar,:,:)-irf(idxlambda,:,:)+...
        (1+omega_y)*(irf(idxY,:,:)-irf(idxhZ,:,:))+irf(idxDelta,:,:);
    % add wb
    irf(idxwb,:,:) = irf(idxw,:,:);
    % add ws
    irf(idxws,:,:) = irf(idxw,:,:);
    if isNatVars
        % add real lending rate (natural)
        irf(idxRrbn,:,:) = irf(idxRrdn,:,:);
        % add cb (natural)
        irf(idxcbn,:,:) = irf(idxhCbar,:,:)-sigma*irf(idxlambdan,:,:);
        % add cs (natural)
        irf(idxcsn,:,:) = irf(idxhCbar,:,:)-sigma*irf(idxlambdan,:,:);
        % add w (natural)
        irf(idxwn,:,:) = irf(idxhmu_w,:,:)-nu*irf(idxhHbar,:,:)-irf(idxlambdan,:,:)+...
            (1+omega_y)*(irf(idxYn,:,:)-irf(idxhZ,:,:));
        % add wb (natural)
        irf(idxwbn,:,:) = irf(idxwn,:,:);
        % add ws (natural)
        irf(idxwsn,:,:) = irf(idxwn,:,:);
    end
    % save the irf
    IRF.(Policy{jPol}) = permute(irf,[1,3,2]);
end

%% ------------------------------------------------------------------------

%% save information
ExerciseName = ['Output_',ExerciseName];
fprintf('\nMAT file: %s\n',ExerciseName)
save(ExerciseName)

%% Elapsed time
disp(' '), vctoc(ttic), disp(' ')
% diary off

%% ------------------------------------------------------------------------
