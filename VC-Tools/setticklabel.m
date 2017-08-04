function [XTick,XTickLabel] = setticklabel(tid,varargin)

% setticklabel
%
% Sets tick and tick labels from list of dates tid.
%
% Usage
%   [XTickLabel,XTick] = setticklabel(tid)
%   Labels will default to first and last year found in tid.
%
%   [XTickLabel,XTick] = setticklabel(...,op)
%   [XTickLabel,XTick] = setticklabel(...,OpName,OpValue)
%   Options abvailable:
%
%     Labels (cell array)
%     List of dates to show in labels. Default: first quarter of first and last 
%     years in tid.
%
%     AnnualTicks (logical) 
%     If set to 1, Ticks set to 1st quarter of each year. If set to 0, unless 
%     ticks is specified on input and is not empty, ticks are set to dates of 
%     labels. Default = 0.
%
%     HideQuarters (logical) 
%     If set to 1 then quarter is hidden from labels. Default: 1.
%
% See also:
% timeidx, findquarter
% 
% ...........................................................................
% 
% Created: April 15, 2017 by Vasco Curdia
% 
% Copyright 2017 by Vasco Curdia

%% options
op.Labels = [];
op.Ticks = [];
op.AnnualTicks = 0;
op.HideQuarters = 1;

op = updateoptions(op,varargin{:});

%% Use only labels present in tid
if isempty(op.Labels)
    q1 = findquarter(tid);
    labels = tid(q1([1,end]));
else
    labels = tid(ismember(tid,op.Labels));
end

%% set tick dates
if isempty(op.Ticks)
    ticks = labels;
    if op.AnnualTicks
        for j=1:4
            idx = findquarter(labels,j);
            tf = all(idx);
            if tf
                q = j;
                break
            end
        end
        if tf
            ticks = tid(findquarter(tid,q));
        end
    end
else
    ticks = tid(ismember(tid,op.Ticks));
end

%% Check labels
if all(ismember(labels,ticks))
    XTick = find(ismember(tid,ticks));
    XTickLabel = ticks;
    XTickLabel(~ismember(XTickLabel,labels)) = {''};
    if op.HideQuarters
        hideq = @(s)s(1:min(4,length(s)));
        XTickLabel = cellfun(hideq,XTickLabel,'UniformOutput',0);
    end
else
    XTick = find(ismember(tid,labels));
    XTickLabel = labels;
end

