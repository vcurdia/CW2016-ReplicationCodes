% plotmgutil
%
% This function plots the marginal utilities of the two types of household
% preferences at steady state.
%
%
% by Vasco Curdia 
% Created: June 23, 2008

%% preamble
clear all
tic
setpath
set(0,'defaultTextInterpreter','latex');

%% Reference values
lambda_b_ss = 2.202679;
lambda_s_ss = 1.809179;
c_b_ss = 0.782089;
c_s_ss = 0.617911;
Cbar_b_ss = 42329.106543;
Cbar_s_ss = 3.174485;
sigma_b = 13.801913;
sigma_s = 2.760383;
yBounds = [0, 5];
xBounds = [0, 2];

%% Evaluate MgUtil
c = xBounds(1):0.01:xBounds(2);
MgUtil_b = (c./Cbar_b_ss).^(-1./sigma_b);
MgUtil_s = (c./Cbar_s_ss).^(-1./sigma_s);

%% Plot

figure

% marginal utility functions
plot(c,MgUtil_b,'-b',c,MgUtil_s,'--r','LineWidth',2)
hold on
xlim(xBounds)
ylim(yBounds)
set(gca,'YTick',yBounds(1):1:yBounds(2))
text(xBounds(2)-0.2,MgUtil_b(end)+0.4,'$u^b_c(c)$','Color','b','FontSize',12)
text(xBounds(2)-0.2,MgUtil_s(end)-0.3,'$u^s_c(c)$','Color','r','FontSize',12)

% steady state levels of consumption
line(c_b_ss*[1;1],[yBounds(1);lambda_b_ss],...
     'Color','b','LineStyle',':','LineWidth',1)
line(c_s_ss*[1;1],[yBounds(1);lambda_s_ss],...
     'Color','r','LineStyle',':','LineWidth',1)
text(c_b_ss+0.025,0.25,'$\bar{c}_b$','Color','b','FontSize',12)
text(c_s_ss-0.075,0.25,'$\bar{c}_s$','Color','r','FontSize',12)

% steady state levels of marginal utility
line([xBounds(1);c_b_ss],lambda_b_ss*[1;1],...
     'Color','b','LineStyle',':','LineWidth',1)
line([xBounds(1);c_s_ss],lambda_s_ss*[1;1],...
     'Color','r','LineStyle',':','LineWidth',1)
text(-0.1,lambda_b_ss+0.03,'$\bar{\lambda}_b$','Color','b','FontSize',14)
text(-0.1,lambda_s_ss-0.1,'$\bar{\lambda}_s$','Color','r','FontSize',14)

%% make eps
print('-dpdf','Fig_01_PlotMgUtil.pdf')

