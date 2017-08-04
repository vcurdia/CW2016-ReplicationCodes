function vcbivhist(x,varargin)

% vcbivhist
%
% Creates the bivariate histogram between of two variables
%
% Created: April 24, 2012 by Vasco Curdia
% Updated: April 25, 2012 by Vasco Curdia
%          Default bins use percentiles to generate bin edges.
%          Option to show marginal means and medians (off by default).
%
% Copyright 2012 by Vasco Curdia

%% Defaults
% xBins = 100*[1,1];
ShowMean = 0;
ShowMedian = 0;
MeanColor = 'c';
MedianColor = 'g';
% LineHeight = 1;
LineWidth = 1;
LineEraseMode = 'normal';
EdgePrc = [0.1,99.9];
nBins = 200;

%% Check options
if ~isempty(varargin)
    nOptions = length(varargin);
    if mod(nOptions,2), error('Incorrect number of optional arguments.'), end
    for jO=1:nOptions/2
        eval(sprintf('%s = varargin{%.0f};',varargin{(jO-1)*2+1},jO*2))
    end
end

%% Run
if ~exist('ParEdges','var')
  xPrc = prctile(x,EdgePrc,1);
  ParEdges = {...
    xPrc(1,1):(xPrc(2,1)-xPrc(1,1))/nBins:xPrc(2,1),...
    xPrc(1,2):(xPrc(2,2)-xPrc(1,2))/nBins:xPrc(2,2),...
    };
end
hist3(x,'Edges',ParEdges);
% hist3(x,xBins);

%     surf(Cj{1},Cj{2},Nj')
set(gcf,'renderer','opengl');
set(get(gca,'child'),'FaceColor','interp','CDataMode','auto');
set(gca,'Color','w');
set(get(gca,'child'),'EdgeColor','interp');
view([0,90])
% %     colormap(vcColorScheme('BaseColors',ColorList(1,:),'LightFactors',1:-0.02:0))
%     colormap(vcColorScheme('BaseColors',[0,0,0],'LightFactors',1:-0.02:0))
%     colormap(vcColorScheme('BaseColors',[0,0,1],'LightFactors',1:-0.02:0))
%     colormap(vcColorScheme('BaseColors',[155,187,89]/256*0.85,'LightFactors',1:-0.02:0))
colormap(hot)
%     colormap(vcColorScheme('BaseColors',ColorList(jc,:),'LightFactors',1:-0.02:0))
%     ColorStack = [ColorStack;vcColorScheme('BaseColors',ColorList(jc,:),'LightFactors',1:-0.02:0)];
%     colormap(ColorStack)
%     cmap = colormap(hot);
%     colormap(cmap(end:-1:1,:))
%     colorbar
%     if exist('ParEdges','var')
%       [Nj,Cj]=hist3(d(j).xd','Edges',ParEdges);
%     else
%       [Nj,Cj]=hist3(d(j).xd',100*[1,1]);
%     end
%     [cc,hh(j)] = contour(Cj{1},Cj{2},Nj,3,'Color',ColorList(j,:),'LineWidth',2);
%     xlim([0,3])
%     ylim([0,10])
axis tight

%% Add mean and median
if ~exist('LineHeight','var')
  xN = hist3(x,'Edges',ParEdges);
  LineHeight = max(xN(:))/2;
end
if ShowMean
  xMean = mean(x,1);
  line(xMean(1)*[1,1],ylim,LineHeight*[1,1],...
    'Color',MeanColor,'LineWidth',LineWidth,'EraseMode',LineEraseMode)
  line(xlim,xMean(2)*[1,1],LineHeight*[1,1],...
    'Color',MeanColor,'LineWidth',LineWidth,'EraseMode',LineEraseMode)
end
if ShowMedian
  xMedian = median(x,1);
  line(xMedian(1)*[1,1],ylim,LineHeight*[1,1],...
    'Color',MedianColor,'LineWidth',LineWidth,'EraseMode',LineEraseMode)
  line(xlim,xMedian(2)*[1,1],LineHeight*[1,1],...
    'Color',MedianColor,'LineWidth',LineWidth,'EraseMode',LineEraseMode)
end

