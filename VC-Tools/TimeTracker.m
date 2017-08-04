classdef (ConstructOnLoad = 1) TimeTracker < handle

    properties
        MasterTimer
        Timers
        ShowOnTimerStop = 1;
    end
    
    methods
        
        function obj = TimeTracker
            obj.MasterTimer = tic;
        end
        
        function start(obj,tn)
            obj.Timers.(tn) = tic;
        end
        
        function stop(obj,tn)
            if ~isfield(obj.Timers,tn)
                error('Timer %s not initiated.',tn)
            end
            if isa(obj.Timers.(tn),'double')
                error('Timer %s was stopped previously.',tn)
            end
            obj.Timers.(tn) = toc(obj.Timers.(tn));
            if obj.ShowOnTimerStop
                obj.showtimers({tn})
            end
        end
        
        function show(obj)
            fprintf('%s\n',vctoc(obj.MasterTimer))
        end
        
        function showtimers(obj,tList)
            if isempty(obj.Timers)
                fprintf('There are no timers to be shown.\n')
                return
            end
            if nargin<2 || ismember('All',tList)
                tList = fieldnames(obj.Timers);
            end
            strLength = int2str(max(strlength(tList)));
            for j=1:length(tList)
                tn = tList{j};
                fprintf(['%',strLength,'s: %s\n'],tn,vctoc(obj.Timers.(tn)))
            end
%             fprintf('\n')
        end
        
    end % methods
    
    methods(Static)
        function obj = loadobj(obj)
            obj.MasterTimer = tic;
        end
    end
    
    
end % classdef
