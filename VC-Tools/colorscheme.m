function cList = colorscheme(varargin)

% colorscheme
%
% ..............................................................................
%
% Created: March 1, 2009 by Vasco Curdia
% 
% Copyright 2009-2017 by Vasco Curdia

%% Generate Color Scheme

%% ------------------------------------------------------------------------

ColorOptions.OldPreferred = [...
         0.15         0.25    0.75; % Blue
         0    0.5000         0; % Green
    0.75         0.15         0.15; % red
    0.2500    0.2500    0.2500; % Grey
%          0    0.5500    0.5500; % Dark Cyan
         0.9,0.6,0.05; %orange
%          0    0.7500    0.7500;
%     0 0 0; % black
%   0.5,0.20,0.7; % purple
%     0.7500         0    0.7500; % Magenta
%     0.7500    0.7500         0;
%     0.2500    0.2500    0.2500;
    ];

ColorOptions.Matlab = [...
         0         0    1;
         0    0.5000         0;
    1         0         0;
    0.2500    0.2500    0.2500; % Grey
         0    0.5500    0.5500; % Dark Cyan
         0.9,0.6,0.05; %orange
%          0    0.7500    0.7500;
%     0 0 0; % black
%   0.5,0.20,0.7; % purple
%     0.7500         0    0.7500; % Magenta
%     0.7500    0.7500         0;
    0.2500    0.2500    0.2500;
    ];

ColorOptions.Brighter = [...
  0.1,0.1,1; % blue
  1,0.1,0.1; % red
%   0,0,0; % black
  0.1,0.8,0.2; % green
  1,0.7,0.05; % orange
  0.1,0.8,0.8; % cyan
  0.5,0.20,0.7; % purple
  ];

ColorOptions.Excel2010 = [[...
    79   129   189;
   192    80    77;
   155   187    89;
   247   150    70;
%    128   100   162;
%     75   172   198;
%    247   150    70;
   ]./256*0.85;
   0.25 0.25 0.25; % add grey
   ];

BaseFactor = 1;
LightFactors = [0,0.25,0.5];
isShowSchemePlot = 0;
isNewFigure = 0;
isLinePlot = 0;
LineWidth = 2;
UseColors = 'Excel2010';

%% Update options
for jO=1:(length(varargin)/2) 
  eval(sprintf('%s = varargin{%.0f};',varargin{(jO-1)*2+1},jO*2))
end

%% preliminary stuff
BaseColors = ColorOptions.(UseColors);
nBaseColors = size(BaseColors,1);
nLightFactors = length(LightFactors);
if ~exist('nColors','var')
    nColors = nBaseColors*nLightFactors;
end

%% Prepare Scheme
BaseColors = BaseFactor*BaseColors;
ColorGap = 1-BaseColors;
for j=1:nLightFactors
    ColorList((j-1)*nBaseColors+(1:nBaseColors),:) = BaseColors+LightFactors(j)*ColorGap;
end
cList = ColorList(1:nColors,:);

%% Plot
if isShowSchemePlot && ~isLinePlot
    if isNewFigure,figure,end
    plot(0:1,0:1,'w')
    hold on
    for j=1:nColors
        fill([j-1,j-1,j,j]/nColors,[0,0.5+j/nColors/2,0.5+j/nColors/2,0],cList(j,:))
    end
    hold off
elseif  isShowSchemePlot && isLinePlot
    if isNewFigure,figure,end
    plot(0:1,0:1,'w')
    hold on
    for j=1:nColors
        plot(0:.1:1,j/(nColors+1)/3+(0:.1:1)*j/(nColors+1),'Color',cList(j,:),'LineWidth',LineWidth)
    end
    hold off
end
