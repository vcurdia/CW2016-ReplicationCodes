function tightfig(h,shape,ax,varargin)

% tightfig
%
% Reduces whitespace around axes in the figure.
%
% h: figure handle
% shape: shape of subplot as cell array
% ax: axes array in same order as subplot
%
% Options:
%   Slack
%   additional slack at edges. default: 0.005
%
% See also
% vcfigure
%
% .........................................................................
%
% Created: April 18, 2017 by Vasco Curdia
% 
% Copyright 2017 by Vasco Curdia


%% options
op.Slack = 0.005;

op = updateoptions(op,varargin{:});

%% existing positions and dimensions
nPlots = length(ax);

axPos = zeros(nPlots,4);
axOutPos = zeros(nPlots,4);
axTI = zeros(nPlots,4);
for jPlot=1:nPlots
    axPos(jPlot,:) = ax(jPlot).Position; 
    axOutPos(jPlot,:) = ax(jPlot).OuterPosition; 
    axTI(jPlot,:) = ax(jPlot).TightInset; 
end

xMinOld = min(axPos(:,1)-axTI(:,1));
xMaxOld = max(axPos(:,1)+axPos(:,3)+axTI(:,3));
wOld = (xMaxOld-xMinOld)/min(nPlots,shape{2});
xPosMinOld = min(axPos(:,1));

yMinOld = min(axPos(:,2)-axTI(:,2));
yMaxOld = max(axPos(:,2)+axPos(:,4)+axTI(:,4));
hOld = (yMaxOld-yMinOld)/min(ceil(nPlots/shape{2}),shape{1});
yPosMinOld = min(axPos(:,2));

%% Detect legend, if below plots
dLeg = 0;
hLeg = findobj(h,'type','legend');
if ~isempty(hLeg)
    legPos = hLeg.Position;
    if legPos(2)<min(axOutPos(:,2))
        dLeg = legPos(4);
    end
end

%% New scale and positions
wNew = (1-2*op.Slack)/shape{2};
hNew = (1-dLeg-2*op.Slack)/shape{1};
[~,xMinIdx] = min(axPos(:,1)-axTI(:,1));
xMin = op.Slack+axTI(xMinIdx,1);
[~,yMinIdx] = min(axPos(:,2)-axTI(:,2));
yMin = op.Slack+dLeg+axTI(yMinIdx,2)+(shape{1}-ceil(nPlots/shape{2}))*hNew;

%% Adjust plot positions
wRatio = wNew/wOld;
hRatio = hNew/hOld;
for jR=1:shape{1}
    for jC=1:shape{2}
        jPlot = (jR-1)*shape{2}+jC;
        if jPlot<=nPlots
            ax(jPlot).Position = [...
                xMin+(axPos(jPlot,1)-xPosMinOld)*wRatio,...
                yMin+(axPos(jPlot,2)-yPosMinOld)*hRatio,...
                axPos(jPlot,3)*wRatio,axPos(jPlot,4)*hRatio];
        end
    end
end
