function IntModelNoFF(varargin)

% IntModelNoFF
%
% Solves for LQ optimal policy in Curdia and Woodford (2008)
% In this version we set financial frictions to zero.
%
% Usage:
%   IntModelNoFF
%   IntModelNoFF(...,OptionName,...)
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
% Created: January 28, 2008
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
ExerciseName = 'NoFF';
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
delta = 0.975;
pi_b = 0.5;
rho_b = 3.2;
Yss = 1;
s_c = 0.7;
% s_bs = 1.2657; % 1.2657
% s_s = s_c/(pi_b*s_bs+1-pi_b);
% s_b = s_bs*s_s;
sigma_bs = isSmSigma+~isSmSigma*5;
psi = 1;
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
    % regular variables
    'Rd', 1.010000;
    'Pi', 1.000000;
    'Y', 0.269529;
    'lambda_b', 1.034759;
    'lambda_s', 1.026946;
    'b', 0.053906;
    'Delta', 1.000000;
    'K', 0.636985;
    'F', 0.636985;
    'b_gy', 0.000000;
    };
if isNatVars
    y = cat(1,y,{...
        % natural variables
        'Rrdn', 1.010000;
        'Yn', 0.269529;
        'lambda_bn', 1.034759;
        'lambda_sn', 1.026946;
        'bn', 0.053906;
        });
end

%% Identify which variables are to be log-linearized
% yLogIdx = ~ismember(y(:,1),'b');
% yLogIdx = true(1,size(y,1));
yLogIdx = ~ismember(y(:,1),{'b_gy'});

%% Exogenous variables
csi = {'xi_i','hZ','hmu_w','htau','hG','hb_g','hHbar','hCbar_b','hCbar_s'};
S = diag([rho_m,rho_csi*ones(1,length(csi)-1)]); 

%% Lagrange multipliers for F constraints
nF = 5+isNatVars*3; for jF=1:nF,FLM(jF,:)={sprintf('FLM%.0f',jF),0};end

%% Lagrange multipliers for G constraints
nG = 4+isNatVars*2; for jG=1:nG,GLM(jG,:)={sprintf('GLM%.0f',jG),0};end

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
Omega_ss = 1;
psi_bs = Omega_ss;
psi_s = psi*(pi_b*psi_bs^(-1/nu)+(1-pi_b))^nu;
psi_b = psi_bs*psi_s;
lambda_s_ss = mu_p*phi*mu_w_ss*Hbar_ss^(-nu)*(Y_ss/Z_ss)^(1+omega_y)/...
    ((1-tau_ss)*Y_ss*(pi_b*psi_b^(-1/nu)+(1-pi_b)*psi_s^(-1/nu))^nu);
lambda_b_ss = lambda_s_ss;
Lambda_ss = lambda_s_ss;
lambdatil_ss = lambda_s_ss;
Lambdatil_ss = lambda_s_ss;
b_ss = rho_b*Y_ss;
sb_ss = (1-delta*Rd_ss)/pi_b/(1-pi_b)*rho_b+...
    (1-delta*Rd_ss)/(1-pi_b)*rho_bg+...
    ((Omega_ss/psi_bs)^(1/nu)-1)/(pi_b*(Omega_ss/psi_bs)^(1/nu)+(1-pi_b))*...
        (1-tau_ss)/mu_p/phi;
s_b = s_c+(1-pi_b)*sb_ss;
s_s = s_c-pi_b*sb_ss;
s_bs = s_b/s_s;
s_g = 1-s_c;
sigma_s = sigmabar/(pi_b*s_b*sigma_bs+(1-pi_b)*s_s);
sigma_b = sigma_bs*sigma_s;
sigma = sigmabar/s_c;
Cbar_b_ss = s_b*Y_ss*lambda_b_ss^sigma_b;
Cbar_s_ss = s_s*Y_ss*lambda_s_ss^sigma_s;
ctil_b_ss = s_b*Y_ss;
ctil_s_ss = s_s*Y_ss;
b_gy_ss = rho_bg*Y_ss;
K_ss = Lambda_ss*mu_p*phi*psi*mu_w_ss/lambdatil_ss*Hbar_ss^(-nu)*(Y_ss/Z_ss)^(1+omega_y)/...
    (1-alpha*beta);
F_ss = Lambda_ss*(1-tau_ss)*Y_ss/(1-alpha*beta);
Rrdn_ss = Rd_ss/Pi_ss;
Yn_ss = Y_ss;
lambda_bn_ss = lambda_b_ss;
lambda_sn_ss = lambda_s_ss;
bn_ss = b_ss;

% a couple more ratios
s_hb = (psi/psi_b*lambda_b_ss/lambdatil_ss)^(1/nu);
s_hs = (psi/psi_s*lambda_s_ss/lambdatil_ss)^(1/nu);
s_hbs = s_hb/s_hs;

% shocks
Z_t = Z_ss*exp(hZ_t);
mu_w_t = mu_w_ss*exp(hmu_w_t);
tau_t = tau_ss+tau_ss*htau_t;
G_t = s_g*Y_ss+Y_ss*hG_t;
b_g_t = rho_bg*Y_ss+Y_ss*hb_g_t;
Hbar_t = Hbar_ss*exp(hHbar_t);
Cbar_b_t = Cbar_b_ss*exp(hCbar_b_t);
Cbar_s_t = Cbar_s_ss*exp(hCbar_s_t);

% auxiliary definitions
Lambda_t = pi_b*lambda_b_t+(1-pi_b)*lambda_s_t;
lambdatil_t = psi*(pi_b*(lambda_b_t/psi_b)^(1/nu)+(1-pi_b)*(lambda_s_t/psi_s)^(1/nu))^nu;
Lambdatil_t = psi^(1/(1+nu))*(pi_b*psi_b^(-1/nu)*lambda_b_t^((1+nu)/nu)+...
    (1-pi_b)*psi_s^(-1/nu)*lambda_s_t^((1+nu)/nu))^(nu/(1+nu));
ctil_b_t = Cbar_b_t*lambda_b_t^(-sigma_b);
ctil_s_t = Cbar_s_t*lambda_s_t^(-sigma_s);
B_t = ctil_b_t-ctil_s_t-((lambda_b_t/psi_b)^(1/nu)-(lambda_s_t/psi_s)^(1/nu))...
    *(lambdatil_t/psi)^(-(1+nu)/nu)*mu_w_t*Hbar_t^(-nu)*(Y_t/Z_t)^(1+omega_y)*Delta_t;

if isNatVars
    % aux defs for natural vars
    Lambdan_t = pi_b*lambda_bn_t+(1-pi_b)*lambda_sn_t;
    lambdatiln_t = psi*(pi_b*(lambda_bn_t/psi_b)^(1/nu)+(1-pi_b)*(lambda_sn_t/psi_s)^(1/nu))^nu;
    Lambdatiln_t = psi^(1/(1+nu))*(pi_b*psi_b^(-1/nu)*lambda_bn_t^((1+nu)/nu)+...
        (1-pi_b)*psi_s^(-1/nu)*lambda_sn_t^((1+nu)/nu))^(nu/(1+nu));
    ctil_bn_t = Cbar_b_t*lambda_bn_t^(-sigma_b);
    ctil_sn_t = Cbar_s_t*lambda_sn_t^(-sigma_s);
    Bn_t = ctil_bn_t-ctil_sn_t-((lambda_bn_t/psi_b)^(1/nu)-(lambda_sn_t/psi_s)^(1/nu))...
        *(lambdatiln_t/psi)^(-(1+nu)/nu)*mu_w_ss*Hbar_t^(-nu)*(Yn_t/Z_t)^(1+omega_y);
end

%% Utility function
U = pi_b*ctil_b_t^(1-1/sigma_b)*Cbar_b_t^(1/sigma_b)/(1-1/sigma_b)...
    +(1-pi_b)*ctil_s_t^(1-1/sigma_s)*Cbar_s_t^(1/sigma_s)/(1-1/sigma_s)...
    -psi/(1+nu)*(lambdatil_t/Lambdatil_t)^(-(1+nu)/nu)*Hbar_t^(-nu)...
    *(Y_t/Z_t)^(1+omega_y)*Delta_t;

%% F constraints
F = vpa([...
	pi_b*(1-pi_b)*B_t-pi_b*b_g_t+delta*(b_tL+pi_b*b_gy_tL)*Rd_tL/Pi_t-b_t;
    pi_b*Cbar_b_t*lambda_b_t^(-sigma_b)+(1-pi_b)*Cbar_s_t*lambda_s_t^(-sigma_s)+G_t-Y_t;
	alpha*Delta_tL*Pi_t^(theta*(1+omega_y))+...
        (1-alpha)*((1-alpha*Pi_t^(theta-1))/(1-alpha))^(theta*(1+omega_y)/(theta-1))-Delta_t;
	(F_t/K_t)^((theta-1)/(1+omega_y*theta))-(1-alpha*Pi_t^(theta-1))/(1-alpha);
	b_g_t-b_gy_t;
    ]);
if isNatVars
    F = vpa([F;
        % natural variables equations
        pi_b*(1-pi_b)*Bn_t-pi_b*b_g_t+delta*(b_tL+pi_b*b_gy_tL)*Rd_tL/Pi_t-bn_t;
        % pi_b*(1-pi_b)*Bn_t-pi_b*b_g_t+delta*(bn_tL+pi_b*b_gy_tL)*Rrdn_tL-bn_t;
        pi_b*Cbar_b_t*lambda_bn_t^(-sigma_b)+(1-pi_b)*Cbar_s_t*lambda_sn_t^(-sigma_s)+G_t-Yn_t;
        1/lambdatiln_t*mu_p*phi*psi*mu_w_ss*Hbar_t^(-nu)*(Yn_t/Z_t)^(1+omega_y)-(1-tau_ss)*Yn_t;
        ]);
end

%% G constraints
G = vpa([...
	Rd_t*beta*...
        ((delta+(1-delta)*pi_b)*lambda_b_tF/Pi_tF+(1-delta)*(1-pi_b)*lambda_s_tF/Pi_tF)-lambda_b_t;
	Rd_t*beta*((1-delta)*pi_b*lambda_b_tF/Pi_tF+...
        (delta+(1-delta)*(1-pi_b))*lambda_s_tF/Pi_tF)-lambda_s_t;
	Lambda_t/lambdatil_t*mu_p*phi*psi*mu_w_t*Hbar_t^(-nu)*(Y_t/Z_t)^(1+omega_y)+...
        alpha*beta*Pi_tF^(theta*(1+omega_y))*K_tF-K_t;
	Lambda_t*(1-tau_t)*Y_t+alpha*beta*Pi_tF^(theta-1)*F_tF-F_t;
    ]);
if isNatVars
    G = vpa([G;
        % natural variables equations
        Rrdn_t*beta*((delta+(1-delta)*pi_b)*lambda_bn_tF+(1-delta)*(1-pi_b)*lambda_sn_tF)-lambda_bn_t;
        Rrdn_t*beta*((1-delta)*pi_b*lambda_bn_tF+(delta+(1-delta)*(1-pi_b))*lambda_sn_tF)-lambda_sn_t;
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
xdisp({'delta','pi_b','beta',...
    'rho_b','s_g','rho_bg','s_c','s_b','s_s','s_bs','s_hb','s_hs','s_hbs',...
    'sigmabar','sigma','sigma_b','sigma_s','sigma_bs',...
    'Cbar_b_ss','Cbar_s_ss','ctil_b_ss','ctil_s_ss',...
    'Omega_ss','psi','psi_b','psi_s','psi_bs',...
    'lambdatil_ss','Lambda_ss','Lambdatil_ss'...
    'phi_pi','phi_y'...
    'Phi','xi','kappa'})
disp(' ')

%% ------------------------------------------------------------------------

%% LQ Optimal Policy
fprintf('\nSolving for LQ Optimal Policy...\n')

%% Run LQ function to get symbolic matrices
[LQmat,LQk_t]=LQ(y_t,csi_t,U,F,G,beta,S,y_ss,FLM_ss,GLM_ss,yLogIdx);

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
g_t = pi_b*s_b*hCbar_b_t+(1-pi_b)*s_s*hCbar_s_t+hG_t;
g_tL = pi_b*s_b*hCbar_b_tL+(1-pi_b)*s_s*hCbar_s_tL+hG_tL;
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
    irf(idxcb,:,:) = irf(idxhCbar_b,:,:)-sigma_b*irf(idxlambda_b,:,:);
    % add cs
    irf(idxcs,:,:) = irf(idxhCbar_s,:,:)-sigma_s*irf(idxlambda_s,:,:);
    % add w
    irf(idxw,:,:) = irf(idxhmu_w,:,:)-nu*irf(idxhHbar,:,:)-...
        (pi_b*(psi/psi_b*lambda_b_ss/lambdatil_ss)^(1/nu)*irf(idxlambda_b,:,:)+...
        (1-pi_b)*(psi/psi_s*lambda_s_ss/lambdatil_ss)^(1/nu)*irf(idxlambda_s,:,:))+...
        (1+omega_y)*(irf(idxY,:,:)-irf(idxhZ,:,:))+irf(idxDelta,:,:);
    % add wb
    irf(idxwb,:,:) = irf(idxw,:,:)+...
        (1-pi_b)/nu*(psi/psi_s*lambda_s_ss/lambdatil_ss)^(1/nu)*...
        (irf(idxlambda_s,:,:)-irf(idxlambda_b,:,:));
    % add ws
    irf(idxws,:,:) = irf(idxw,:,:)-...
        pi_b/nu*(psi/psi_b*lambda_b_ss/lambdatil_ss)^(1/nu)*...
        (irf(idxlambda_s,:,:)-irf(idxlambda_b,:,:));
    if isNatVars
        % add real lending rate (natural)
        irf(idxRrbn,:,:) = irf(idxRrdn,:,:);
        % add cb (natural)
        irf(idxcbn,:,:) = irf(idxhCbar_b,:,:)-sigma_b*irf(idxlambda_bn,:,:);
        % add cs (natural)
        irf(idxcsn,:,:) = irf(idxhCbar_s,:,:)-sigma_s*irf(idxlambda_sn,:,:);
        % add w (natural)
        irf(idxwn,:,:) = irf(idxhmu_w,:,:)-nu*irf(idxhHbar,:,:)-...
            (pi_b*(psi/psi_b*lambda_b_ss/lambdatil_ss)^(1/nu)*irf(idxlambda_bn,:,:)+...
            (1-pi_b)*(psi/psi_s*lambda_s_ss/lambdatil_ss)^(1/nu)*irf(idxlambda_sn,:,:))+...
            (1+omega_y)*(irf(idxYn,:,:)-irf(idxhZ,:,:));
        % add wb (natural)
        irf(idxwbn,:,:) = irf(idxwn,:,:)+...
            (1-pi_b)/nu*(psi/psi_s*lambda_s_ss/lambdatil_ss)^(1/nu)*...
            (irf(idxlambda_sn,:,:)-irf(idxlambda_bn,:,:));
        % add ws (natural)
        irf(idxwsn,:,:) = irf(idxwn,:,:)-...
            pi_b/nu*(psi/psi_b*lambda_b_ss/lambdatil_ss)^(1/nu)*...
            (irf(idxlambda_sn,:,:)-irf(idxlambda_bn,:,:));
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
