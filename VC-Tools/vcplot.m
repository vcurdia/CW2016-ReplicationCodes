function h = vcplot(y,varargin)

% vcplot
%
% Plots several lines with optional light effects.
% 
% Syntax:
%   vcplot(y)
%   vcplot(x,y)
%   vcplot(...,OptionsStructure,...)
%   vcplot(...,'PropertyName',PropertyValue,...)
%   h = vcplot(...)
%
% If no arguments are specified then an example is shown.
%
% See also:
% vcfigure, vcplotdistbands, colorscheme
%
% .........................................................................
%
% Created: April 24, 2013
%
% Copyright 2013-2017 by Vasco Rafael da Silva Curdia

%% ---------------------------------------------------------------------

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
elseif nargin==1
    x = 1:size(y,2);
else
    if isnumeric(varargin{1})
        x = y;
        y = varargin{1};
        varargin(1) = [];
    end
end

%% Dimensions
[ny,nx] = size(y);

%% Default options
op.Color = colorscheme;
op.LineStyle = {'-','--','-.',':',':',':'};
op.LineMarker = {'none','none','none','o','s','^'};
% op.LineStyle = {'-','--','-','-','none'};
% op.LineMarker = {'none','none','o','s','^'};
% op.LineMarkerSize = [2,2,2,2,2,2];
% op.LineWidth = [1.5,1.5,1.5,0.5,0.5,0.5];
op.LineMarkerSize = [1,1,1,2,2,2];
op.LineWidth = [1,1,1,0.5,0.5,0.5];
op.ShowZeroLine = 1;
op.ZeroLineColor = 'k';
op.ZeroLineStyle = '-';
op.ZeroLineWidth = 0.5;
op.ShowLegend = (ny>1);
op.LegendString = {};
op.LegendLocation = 'Best';
op.LegendOrientation = 'vertical'; %'vertical','horizontal'
op.LegendItems = [];

%% Update options
op = updateoptions(op,varargin{:});


%% Check x
if ~exist('x','var')
  x = 1:nx;
end

%% Some values
CurrentHold = ishold;
% xx = x(1):1/op.MarkerTicks:x(nx);
xx = x;
nxx = length(xx);

%% Plot zero line
if op.ShowZeroLine
    h.ZeroLine = plot(xx,zeros(1,nxx),op.ZeroLineStyle,...
                      'Color',op.ZeroLineColor,'LineWidth',op.ZeroLineWidth);
    h.ZeroLine.Annotation.LegendInformation.IconDisplayStyle = 'off';
    if ~ishold,hold on,end
else
    h.ZeroLine = [];
end

%% Plot Line
for j=1:ny
    h.Lines(j) = plot(x,y(j,:),...
                      'LineStyle',op.LineStyle{j},...
                      'Color',op.Color(j,:),...
                      'MarkerFaceColor',op.Color(j,:),...
                      'Marker',op.LineMarker{j},...
                      'MarkerSize',op.LineMarkerSize(j),...
                      'LineWidth',op.LineWidth(j));
    if ~ishold,hold on,end
end

%% Show Legend
if op.ShowLegend
    if isempty(op.LegendString)
        for j=1:ny
            op.LegendString{j} = sprintf('Line %.0f',j);
        end
    end
    h.LegendItems = h.Lines;
    if ~isempty(op.LegendItems)
        h.LegendItems = [h.LegendItems,op.LegendItems];
    end
    nLegendItems = length(h.LegendItems);
    h.Legend = legend(h.LegendItems,op.LegendString,...
                      'Location',op.LegendLocation,...
                      'Orientation',op.LegendOrientation);
else
    h.LegendItems = [];
    h.Legend = [];
end

%% Restore hold status
if ~CurrentHold
    hold off
end

%% end example
if nargin==0
    hold off
    axis tight
end

%% Exit
h.Options = op;

if nargout==0
    clear h
end

%% ------------------------------------------------------------------------
