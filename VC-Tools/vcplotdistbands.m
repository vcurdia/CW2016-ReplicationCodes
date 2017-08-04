function h = vcplotdistbands(y,varargin)

% vcplotdistbands
%
% Plots median and percentile bands for matrix X.
%
% Usage: 
%   vcplotdistbands(y)
%   vcplotdistbands(x,y)
%   vcplotdistbands(...,OptionsStructure,...)
%   vcplotdistbands(...,'PropertyName',PropertyValue,...)
%   h = vcplotdistbands(...)
%
% If no arguments are specified then an example is shown.
%
% Required input argument:
%
%   y
%   Matrix containing data. Percentiles computed along first dimension.
%
% Options:
%
%   Bands2Show
%   Percent intervals to be shown in the plots, centered around median.
%   Default: [50,60,70,80,90]
%
%   MedianColor
%   Color of median line.
%   Default: [0,0,0.7]
%
%   ShadeColor
%   Base color for bands.
%   Default: [0.2,0.6,0.5]
%
%   LineWidth
%   Width of median line.
%   Default: 1.5
%
%   isZeroLine
%   If 1 plots the zero line. If 0 it does not plot a zero line.
%   Default: 1
%
%   ZeroLineColor
%   Color for the zero line.
%   Default: 'k'
%
%   ZeroLineStyle
%   Style for the zero line.
%   Default: ':'
%
%   tid
%   x axis values.
%   Default: 1:T
%
% See also:
% vcplot, vcfigure, colorscheme
%
% ..............................................................................
%
% Created: October 30, 2008 by Vasco Curdia
% 
% Copyright 2008-2017 by Vasco Curdia

%% ------------------------------------------------------------------------

%% Check Inputs
if nargin==0
% Example
    x = 1:25;
    y = rand(1000,length(x));
    op.AltData = [...
        -0.1+x/25;
        0.4+x/25;
%     1.3-x/25;
                ];
    op.LegendLocation = 'SE';
    varargin{1} = op;
else
    if isnumeric(varargin{1})
        x = y;
        y = varargin{1};
        varargin(1) = [];
    end
end

%% default options
op.Bands2Show = [50,70,90]; %[50,60,70,80,90]
op.AltData = [];
op.ShowLegend = 0;
op.LegendLocation = 'Best';
op.LegendString = {};
op.LegendOrientation = 'vertical'; %'vertical','horizontal'
op.LegendWithBands = 0;
op.LineColor = colorscheme;
% op.LineColor = op.LineColor([3,1,2,4:end],:); 
  % Assumes that 3rd line color is red and first is blue
% op.ShadeColor = [0.2,0.6,0.5];
% op.ShadeColor = [0.45,0.45,0.5]; 
% op.ShadeColor = [0.585,0.585,0.65]; 
op.ShadeColor = [0.72,0.77,0.82]*0.95;
% op.ShadeColor = [0.15,0.25,0.75];
op.ShadeColorBrightness = 0.9;
% op.ShadeFactors = [0.1,0.65]% shade factors at 50 and 90%
% op.ShadeFactors = [0.2,0.7]; % shade factors at 50 and 90%
op.ShadeFactors = [0.1,0.65]; % shade factors at 50 and 90%
% op.ShadeFactors = [0.1,0.7]; % shade factors at 50 and 90%
% MedianColor = [0,0,0.7];
% ShadeColor = [0.2,0.6,0.5];
op.FaceAlpha = 1;

%% Update options
op = updateoptions(op,varargin{:});

%% Check options
% op.LineColor = [0.3,0.75,0.3;op.LineColor];
% if ~isfield(op,'ShadeColor')
%     if ~isfield(op,'ShadeColorBase')
% %   op.ShadeColorBase = [0.70,0.70,0.7]; % grey
% %   op.ShadeColorBase = [0.70,0.80,0.85]; % light blue 
% %   op.ShadeColorBase = [0.60,0.70,0.75]; % grey/Blue
%         op.ShadeColorBase = [0.15,0.25,0.75]; % Blue
%     end
%     if ~isfield(op,'ShadeCLWeight')
%         op.ShadeCLWeight = 0.0;
%     end
% %   op.ShadeColor = [0.65,0.75,0.8]*0.95+0.05*op.LineColor(1,:);
%     op.ShadeColor = ...
%         op.ShadeColorBrightness*op.ShadeColorBase*(1-op.ShadeCLWeight)+...
%         op.ShadeCLWeight*op.LineColor(1,:);
% end

%% Check x
nx = size(y,2);
if ~exist('x','var')
    x = 1:nx;
end

%% Plot bands
op.Bands2Show = sort(op.Bands2Show,'descend');
nBands = length(op.Bands2Show);
InitHold = ishold;
BandsData = zeros(nBands*2,nx);
BandColorSlope = [-1,1]*op.ShadeFactors'/([1,-1]*op.Bands2Show([1,nBands])');
BandColorCt = op.ShadeFactors(1)-BandColorSlope*op.Bands2Show(nBands);
for jB=1:nBands
    Band = op.Bands2Show(jB);
    BandPath = prctile(y,50+Band/2*[-1,+1]);
    BandColor = op.ShadeColor+...
        (1-op.ShadeColor)*(BandColorCt+BandColorSlope*Band);
    h.Bands(jB) = fill([x,x(end:-1:1)],[BandPath(1,:),BandPath(2,end:-1:1)],...
                       BandColor,'EdgeColor',BandColor);
    hold on
    BandsData((jB-1)*2+[1,2],:) = BandPath;
end
h.BandsData = BandsData;

%% Plot lines
YData = [prctile(y,50);op.AltData];
ny = size(YData,1);
if op.ShowLegend
    if isempty(op.LegendString)
        op.LegendString{1} = 'Median';
        for j=2:ny
            op.LegendString{j} = sprintf('Alt %.0f',j-1);
        end
        for j=1:nBands*op.LegendWithBands
            op.LegendString{ny+j} = sprintf('%.0f%%',op.Bands2Show(j));
        end
    end
    if op.LegendWithBands
        op.LegendItems = h.Bands;
    end
end
hp = vcplot(x,YData,op);
h.Lines = hp.Lines;
h.ZeroLine = hp.ZeroLine;
h.Legend = hp.Legend;
h.LegendItems = hp.LegendItems;
h.LineOptions = hp.Options;
alpha(gca,op.FaceAlpha)

%% some more stuff
if ~InitHold
    hold off
end
xlim([x(1),x(end)])


%% Exit
h.XData = x;
h.YData = YData;
h.Options = op;

%% ------------------------------------------------------------------------
