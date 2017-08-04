function varargout = vcplotshadow(x,y,varargin)

% vcplotshadow
%
% Plots several lines with a shadow and a light effect.
% 
% If no arguments are specified then an example is shown.
%
% .........................................................................
%
% Created: April 24, 2013
%
% Copyright (2013) by Vasco Rafael da Silva Curdia

%% ------------------------------------------------------------------------

%% Check main input
if nargin==0
  x = 1:10;
  y(1,:) = -2+x;
  y(2,:) = 2+log(x);
  y(3,:) = 1+2*log(x);
  y(4,:) = 5-2*log(x);
  y(5,:) = 7-2*log(x);
  y(6,:) = 1-2*log(x);
  % y(2,1:1) = NaN;
end

%% Default options
LineColor = vcColorScheme;
LineStyle = {'-','--','-.','none','none','none'};
LineMarker = {'none','none','none','o','s','^'};
LineMarkerSize = [2,2,2,2,2,2];
LineWidth = 2;
ShowShadow = 0;
ShowLight = 1;
UseLightAngle = 0;
ShadowColor = 0.9*[1,1,1];
Brightness = 0.3;
DistFactor = 0.001;
LegendDistFactor = 0.0005;
ShowZeroLine = 1;
ShowLegend = (size(y,1)>1);
ZeroLineColor = 'k';
ZeroLineStyle = ':';
ZeroLineWidth = 0.5;
LegendLocation = 'Best';
LegendOrientation = 'vertical'; %'vertical','horizontal'

%% Update options
if ~isempty(varargin)
  nOptions = length(varargin);
  if nOptions==1 && isstruct(varargin{1})
    Options = fieldnames(varargin{1});
    for jO=1:length(Options)
      eval(sprintf('%1$s = varargin{1}.%1$s;',Options{jO}))
    end
  elseif mod(nOptions,2)
    error('Incorrect number of optional arguments.')
  else
    for jO=1:nOptions/2
      eval(sprintf('%s = varargin{%.0f};',varargin{(jO-1)*2+1},jO*2))
    end
  end
end

%% Some values
ny = size(y,1);
xDist = (max(x)-min(x))*DistFactor;
yDist = (max(y(:))-min(y(:)))*DistFactor;
CurrentHold = ishold;
h.lines = zeros(4,ny);
if ShowLegend && ~exist('LegendString','var')
  for j=1:ny
    LegendString{j} = sprintf('Line %.0f',j);
  end
end

%% Plot Shadow
if ShowShadow
  for j=1:ny
    h.lines(1+2*ShowLight+1,j) = plot(x+2*xDist,y(j,:)-2*yDist,...
      'LineStyle',LineStyle{j},...
      'Color',ShadowColor,'MarkerFaceColor',ShadowColor,...
      'Marker',LineMarker{j},'MarkerSize',LineMarkerSize(j),...
      'LineWidth',LineWidth);
    if ~ishold,hold on,end
  end
end

%% Plot zero line
if ShowZeroLine
  h.zeroline = plot(x,zeros(size(x)),ZeroLineStyle,'Color',ZeroLineColor,...
    'LineWidth',ZeroLineWidth);
  set(get(get(h.zeroline,'Annotation'),'LegendInformation'),...
      'IconDisplayStyle','off');
  if ~ishold,hold on,end
else
  h.zeroline = [];
end

%% Plot Line with light effect
for j=1:ny
  h.lines(1,j) = plot(x,y(j,:),...
    'LineStyle',LineStyle{j},...
    'Color',LineColor(j,:),'MarkerFaceColor',LineColor(j,:),...
    'Marker',LineMarker{j},'MarkerSize',LineMarkerSize(j),...
    'LineWidth',LineWidth);
  if ~ishold,hold on,end
  if ShowLight
    h.lines(2,j) = plot(x-UseLightAngle*xDist,y(j,:)+UseLightAngle*yDist,...
      'LineStyle',LineStyle{j},...
      'Color',LineColor(j,:)+(1-LineColor(j,:))*Brightness/2,...
      'MarkerFaceColor',LineColor(j,:)+(1-LineColor(j,:))*Brightness/2,...
      'Marker',LineMarker{j},'MarkerSize',LineMarkerSize(j)*2/3,...
      'LineWidth',LineWidth*2/3);
    h.lines(3,j) = plot(x-UseLightAngle*2*xDist,y(j,:)+UseLightAngle*2*yDist,...
      'LineStyle',LineStyle{j},...
      'Color',LineColor(j,:)+(1-LineColor(j,:))*Brightness,...
      'MarkerFaceColor',LineColor(j,:)+(1-LineColor(j,:))*Brightness,...
      'Marker',LineMarker{j},'MarkerSize',LineMarkerSize(j)/3,...
      'LineWidth',LineWidth/3);
  end
end

%% Show Legend
if ShowLegend
  [h.legend,h.legendobj,~,~] = legend(LegendString,...
    'Location',LegendLocation,'Orientation',LegendOrientation);
  hlegpos = get(h.legend,'Position');
  xDist = LegendDistFactor/hlegpos(3);
  yDist = LegendDistFactor/hlegpos(4);
  % redo lines
  set(h.legendobj(ny+1:end),'Visible','off')
  if ShowShadow
    for j=1:ny
      line('XData',get(h.legendobj(ny+2*j-1),'XData')+2*xDist,...
        'YData',get(h.legendobj(ny+2*j-1),'YData')-2*yDist,...
        'Parent',h.legend,'LineStyle',LineStyle{j},...
        'Color',ShadowColor,'MarkerFaceColor',ShadowColor,...
        'Marker',LineMarker{j},'MarkerSize',LineMarkerSize(j),...
        'LineWidth',LineWidth);
    end
  end
  for j=1:ny
    line('XData',get(h.legendobj(ny+2*j-1),'XData'),...
      'YData',get(h.legendobj(ny+2*j-1),'YData'),...
      'Parent',h.legend,'LineStyle',LineStyle{j},...
      'Color',LineColor(j,:),'MarkerFaceColor',LineColor(j,:),...
      'Marker',LineMarker{j},'MarkerSize',LineMarkerSize(j),...
      'LineWidth',LineWidth);
    if ShowLight
      line('XData',get(h.legendobj(ny+2*j-1),'XData')-UseLightAngle*xDist,...
        'YData',get(h.legendobj(ny+2*j-1),'YData')+UseLightAngle*yDist,...
        'Parent',h.legend,'LineStyle',LineStyle{j},...
        'Color',LineColor(j,:)+(1-LineColor(j,:))*Brightness/2,...
        'MarkerFaceColor',LineColor(j,:)+(1-LineColor(j,:))*Brightness/2,...
        'Marker',LineMarker{j},'MarkerSize',LineMarkerSize(j)*2/3,...
        'LineWidth',LineWidth*2/3);
      line('XData',get(h.legendobj(ny+2*j-1),'XData')-2*UseLightAngle*xDist,...
        'YData',get(h.legendobj(ny+2*j-1),'YData')+2*UseLightAngle*yDist,...
        'Parent',h.legend,'LineStyle',LineStyle{j},...
        'Color',LineColor(j,:)+(1-LineColor(j,:))*Brightness,...
        'MarkerFaceColor',LineColor(j,:)+(1-LineColor(j,:))*Brightness,...
        'Marker',LineMarker{j},'MarkerSize',LineMarkerSize(j)/3,...
        'LineWidth',LineWidth/3);
    end
  end
else
  h.legend = [];
  h.legendobj =[];
end

%% Restore hold status
if ~CurrentHold
  hold off
end

if nargout>0
  varargout{1} = h;
end

%% end example
if nargin==0
  hold off
  axis tight
  PlotName = 'TestPlot.eps';
  print('-depsc2',PlotName)
  eval(['!epstopdf ',PlotName])
end

%% ------------------------------------------------------------------------
