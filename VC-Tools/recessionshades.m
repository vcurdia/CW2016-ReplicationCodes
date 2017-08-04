function h = recessionshades(Rec,varargin)

% recessionshades
%
% Adds recession shades to a plot. 
%
% Note: This function should be used when the plot is completely finished and no 
% further axis resizing will take place.
%
% Usage:
%
%   recessionshades(Rec)
%   Where Rec is a series with 1 at the positions of the recessions, and 0
%   otherwise.
%
%   recessionshades(...,'TimeIdx',TimeIdx,...)
%   Sets the values time index in the plot. Default: 1,2,3,...
%
%   recessionshades(...,'EdgeColor',value,...)
%   Sets the color of the edge of the shades to value. Default: 'none'
%
%   recessionshades(...,'FaceColor',value,...)
%   Sets the color of the face of the shades to value. Default: .85*[1,1,1]
%
%   recessionshades(...,'LayerPos',value,...)
%   If value is set to 'top' then the axis lines show up on top of the
%   figure. If instead value is set to 'bottom' then the axis lines show up
%   below everything in the figure, including the shades. Default: 'top'
%
% ..............................................................................
% 
% Created: April 26, 2011
% 
% Copyright 2011-2017 by Vasco Curdia

%% -----------------------------------------------------------------------------

%% Defaults
op.EdgeColor = 'none';
op.FaceColor = .85*[1,1,1];
op.LayerPos = 'top';
h.TimeIdx = 1:length(Rec);
op.ShowBaseline = 'off';

%% Check options
op = updateoptions(op,varargin{:});

%% apply shades
yBounds = get(gca,'YLim');
xBounds = get(gca,'XLim');
hold on
h.Shade(1) = bar(h.TimeIdx,yBounds(2)*Rec,1,...
                 'EdgeColor',op.EdgeColor,...
                 'FaceColor',op.FaceColor,...
                 'ShowBaseline',op.ShowBaseline);
h.Shade(2) = bar(h.TimeIdx,yBounds(1)*Rec,1,...
                 'EdgeColor',op.EdgeColor,...
                 'FaceColor',op.FaceColor,...
                 'ShowBaseline',op.ShowBaseline);
uistack(h.Shade(1),'bottom')
uistack(h.Shade(2),'bottom')
hold off
ylim(yBounds)
xlim(xBounds)
set(gca,'Layer',op.LayerPos)

h.Options = op;

%% -----------------------------------------------------------------------------

