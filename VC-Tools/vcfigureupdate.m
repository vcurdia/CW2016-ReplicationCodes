function h = vcfigureupdate(h,NewData,varargin)

% vcfigureupdate
%
% Resuses existing figure. Assumes that the number of subplots and their 
% subtitles are the same.
%
% See also:
% vcfigure, vcplot, vcplotdistbands, vccolorlist, vcRecessionShades
%
% .........................................................................
%
% Created: January 4, 2017
%
% Copyright 2017 by Vasco Curdia

%% ------------------------------------------------------------------------

%% Options
op = h.Options;

%% Preparations
h.YData = NewData;

nPlots = size(NewData,3);
nAltData = size(op.AltData,1);

if h.Options.PlotBands
    Bands = sort(h.Plot(1).Options.Bands2Show,'descend');
    nBands = length(Bands);
    prc = nan(1,1+2*nBands);
    prc(1) = 50;
    for jB=1:nBands
        prc(1+(jB-1)*2+[1,2]) = 50+Bands(jB)/2*[-1,1];
    end
end


%% Update figure
for jP=1:nPlots
    if op.PlotBands
        BandsData = prctile(NewData(:,:,jP),prc,1);
        h.Plot(jP).BandsData = BandsData(2:end,:);
        h.Plot(jP).YData = [BandsData(1,:);op.AltData];
        for jA=1:(1+nAltData)
            h.Plot(jP).Lines(jA).YData = h.Plot(jP).YData(jA,:);
        end
        for jB=1:nBands
            h.Plot(jP).Bands(jB).Vertices(:,2) = ...
                [BandsData(1+(jB-1)*2+1,:),...
                 BandsData(1+(jB-1)*2+2,end:-1:1)]';
        end
        YLim = [min(BandsData(:)),max(BandsData(:))];
    else
        NewDataj = squeeze(NewData(:,:,jP));
%         h.Plot(jP).Lines(1).YData = NewDataj;
        for jL=1:size(NewData,1)
            h.Plot(jP).Lines(jL).YData = NewDataj(jL,:);
        end
        YLim = [min(NewDataj(:)),max(NewDataj(:))];
    end
    yMinScale = op.YMinScale(1+(jP-1)*(length(op.YMinScale)>1));
    YLim = [min([1+op.YSlack,-op.YSlack]*YLim',-yMinScale),...
            max([-op.YSlack,1+op.YSlack]*YLim',yMinScale)];
    if YLim(2)==YLim(1)
        YLim = YLim+[-0.01,+0.01];
    end
    h.SubPlot(jP).YLim = YLim;
    if op.ShowRecessionShades
        h.RecessionShades(jP).Shade(1).YData = ...
            YLim(1)*h.Options.RecessionShades;
        h.RecessionShades(jP).Shade(2).YData = ...
            YLim(2)*h.Options.RecessionShades;
    end
    if op.TightFig
        tightfig(h.Figure,op.Shape,h.SubPlot,op.TightFigOptions)
    end
end


%% ------------------------------------------------------------------------
