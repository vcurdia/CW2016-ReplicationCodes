function h = vcfigure(y,varargin)

% vcfigure
%
% Creates figure with multiple plots (or just one). Plots can have multiple
% lines in each, but all need to have the same number of lines.
%
% Syntax:
%
%   vcfigure(y)
%   vcfigure(x,y)
%   vcfigure(...,OptionsStructure)
%   vcfigure(...,'PropertyName',PropertyValue,...)
%   h = vcfigure(...)
%
% If no arguments are specified then an example is shown.
%
% See also:
% vcplot, vcplotdistbands, vcColorScheme, vcRecessionShades
%
% .........................................................................
%
% Created: October 11, 2013 by Vasco Curdia
%
% Copyright 2013-2017 by Vasco Curdia


%% Options

%% Check Inputs
if nargin==0
% Example with bands, including alternatives
    x = 1:25;
    y = rand(1000,length(x),9);
    op.AltData = repmat([-0.1+x/25;
                        0.4+x/25;
                        1.3-x/25],1,1,9);
    op.PlotBands = 1;
    op.Bands2Show = [50,70,90];
    op.TitleList = {'x1','x2','x3','x4','x5','x6','x7','x8','x9'};
    varargin{1} = op;
% % Example with bands, including alternatives, 3x2
%     x = 1:25;
%     y = rand(1000,length(x),6);
%     op.AltData = repmat([-0.1+x/25;
%                         0.4+x/25;
%                         1.3-x/25],1,1,9);
%     op.PlotBands = 1;
%     op.Bands2Show = [50,70,90];
%     op.Shape = {3,2};
%     op.TitleList = {'x','x','x','x','x','x','x','x','x'};
%     varargin{1} = op;
% % Example with lines
%     x = 1:10;
%     y(1,:) = -2+x;
%     y(2,:) = 2+log(x);
%     y(3,:) = 1+2*log(x);
%     y(4,:) = 5-2*log(x);
%     y(5,:) = 7-2*log(x);
%     y = repmat(y,[1,1,9]);
end
if ~isempty(varargin) && isnumeric(varargin{1})
    x = y;
    y = varargin{1};
    varargin(1) = [];
end

%% Some dimensions
[ny,nx,nPlots] = size(y);

%% Default Options
op.NewFig = 1;
op.Visible = 'on';
op.Shape = {1,1};
op.PlotBands = 0;
op.Plot = struct;
op.Plot.LegendLocation = 'SouthOutside'; % 'SouthOutside','EmptySlot'
op.Plot.LegendOrientation = [];
op.Plot.ShowLegend = [];
op.AxisTight = 1;
op.YMinScale = 0; % Set to 0 if no min slack needed, but want 0 to show. 
                    % Set to NaN in order to not show zero.
op.YSlack = 0.05;
op.AltData = [];
op.RecessionShades = [];
op.RecessionShadesOptions = {};
op.TitleList = {};
op.TitleFontSize = [];
op.XLim = [];
op.XTick = [];
op.XTickLabel = [];
op.YLim = [];
op.YGrid = [];
op.AxesFontSize = [];
op.SupTitle = [];
op.SupTitleFontSize = [];
% op.PaperSize = [6.5, 6.5];
op.PaperSize = [];
op.PaperPosition = [];
op.TightFig = 1;
op.TightFigOptions = struct;

op = updateoptions(op,varargin{:});

%% Check options
if ~isempty(op.PaperSize) && isempty(op.PaperPosition)
    op.PaperPosition = [0, 0, op.PaperSize(1), op.PaperSize(2)];
end
if op.PlotBands
    if ~isempty(op.AltData)
        ny = 1+size(op.AltData,1);
    else
        ny = 1;
    end
end
if isempty(op.Plot.ShowLegend)
%   op.Plot.ShowLegend = (ny>1) && (~o.PlotBands);
    op.Plot.ShowLegend = (ny>1);
end
if isempty(op.Plot.LegendOrientation)
    if nPlots>1 && strcmp(op.Plot.LegendLocation,'SouthOutside')
        op.Plot.LegendOrientation = 'horizontal';
    else
        op.Plot.LegendOrientation = 'vertical';
    end
end
% if ~(nPlots>1 && (strcmp(op.Plot.LegendLocation,'SouthOutside') || ...
%     strcmp(op.Plot.LegendLocation,'EmptySlot')))
%     op.Plot.LegendLocation = op.Plot.LegendLocation;
% end
% if isempty(op.AxesFontSize)
%     if nPlots == 1
%         op.AxesFontSize = 10;
%     elseif nPlots<=4
%         op.AxesFontSize = 8;
%     else
%         op.AxesFontSize = 6;
%     end        
% end


%% Create figure

%% Initiate figure
if op.NewFig
    h.Figure = figure('Visible',op.Visible);
else
    h.Figure = gcf;
    h.Figure.Visible = op.Visible;
    clf
end

%% Check x
if ~exist('x','var')
    x = 1:nx;
end

%% Check Figure Shape
if isempty(op.Shape), op.Shape = {1,1}; end
jDim = 2;
while prod([op.Shape{:}])<nPlots
    op.Shape{jDim} = op.Shape{jDim}+1;
    jDim = ~(jDim-1)+1;
end

%% Check Recession shades
op.ShowRecessionShades = ~isempty(op.RecessionShades);

%% Plot data
for jPlot=1:nPlots
    h.SubPlot(jPlot) = subplot(op.Shape{:},jPlot);
    yj = y(:,:,jPlot);
    opj = op.Plot;
    if ~isempty(op.AltData)
        opj.AltData = op.AltData(:,:,jPlot);
    end
    if op.Plot.ShowLegend && jPlot==nPlots
        opj.ShowLegend = 1;
        if strcmp(op.Plot.LegendLocation,'EmptySlot')
            opj.LegendLocation = 'South';
        end
        if strcmp(op.Plot.LegendLocation,'SouthOutside') && nPlots>1
            opj.LegendLocation = 'South';
        end
    else
        opj.ShowLegend = 0;
    end
    if op.PlotBands
        h.Plot(jPlot) = vcplotdistbands(x,yj,opj);
    else
        h.Plot(jPlot) = vcplot(x,yj,opj);
    end
    if ~isempty(op.TitleList)
        hh = title(op.TitleList{jPlot});
        if ~isempty(op.TitleFontSize)
            hh.FontSize = op.TitleFontSize;
        end
    end
    ax = gca;
    if op.AxisTight, axis tight, end
    if ~isempty(op.XLim)
        ax.XLim = op.XLim;
    end
    if ~isempty(op.XTick)
        ax.XTick = op.XTick;
    end
    if ~isempty(op.XTickLabel)
        ax.XTickLabel = op.XTickLabel;
    end
    if ~isempty(op.YLim)
        if size(op.YLim,3)==1
            yBounds = op.YLim(1,:);
        else
            yBounds = op.YLim(1,:,jPlot);
        end
    else
        yMinScale = op.YMinScale(1+(jPlot-1)*(length(op.YMinScale)>1));
        yBounds = ylim;
        yBounds = [min([1+op.YSlack,-op.YSlack]*yBounds',-yMinScale),...
                   max([-op.YSlack,1+op.YSlack]*yBounds',yMinScale)];
    end
    ylim(yBounds)
    if ~isempty(op.YGrid)
        ax.YGrid = op.YGrid;
    end
    if ~isempty(op.AxesFontSize)
        ax.FontSize = op.AxesFontSize;
    end
    if op.ShowRecessionShades
        h.RecessionShades(jPlot) = ...
            recessionshades(op.RecessionShades,op.RecessionShadesOptions{:});
    end
%     ax = gca;
%     outerpos = ax.OuterPosition;
%     ti = ax.TightInset; 
%     left = outerpos(1) + ti(1);
%     bottom = outerpos(2) + ti(2);
%     ax_width = outerpos(3) - ti(1) - ti(3);
%     ax_height = outerpos(4) - ti(2) - ti(4);
%     ax.Position = [left bottom ax_width ax_height];
end

if ~op.ShowRecessionShades
    h.RecessionShades = [];
end
    
if op.Plot.ShowLegend && nPlots>1
    h.Legend = h.Plot(nPlots).Legend;
    legPos = h.Legend.Position;
    if strcmp(op.Plot.LegendLocation,'SouthOutside')
        legPos(1) = 0.5 - legPos(3)/2;
        legPos(2) = 0;
    elseif strcmp(op.Plot.LegendLocation,'EmptySlot')
        if nPlots<prod([op.Shape{:}])
            xIdx = (max(1,op.Shape{1}-1))*op.Shape{2};
            xR = h.SubPlot(xIdx).Position;
            xL = h.SubPlot(min(nPlots,xIdx+1)).Position;
            legPos(1) = xR(1)+(xR(3)-legPos(3))/2;
            legPos(2) = xL(2)+(xL(4)-legPos(4))/2;
        end
    end
    h.Legend.Position = legPos;
else
    h.Legend = [];
end
if ~isempty(op.SupTitle)
    h.SupTitle = suptitle(op.SupTitle);
    if ~isempty(op.SupTitleFontSize)
        h.SupTitle.FontSize = op.SupTitleFontSize;
    end
end

%% Set Fig paper size for pdf printing
if ~isempty(op.PaperSize)
    h.Figure.PaperSize = op.PaperSize;
end
if ~isempty(op.PaperPosition)
    h.Figure.PaperPosition = op.PaperPosition;
end

%% eliminate needless white space
if op.TightFig
    tightfig(h.Figure,op.Shape,h.SubPlot,op.TightFigOptions)
end


%% Exit
h.XData = x;
h.YData = y;
h.Options = op;

if nargout==0
    clear h
end

