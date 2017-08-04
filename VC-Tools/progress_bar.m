function progress_bar(action,progress,title)
% PROGRESS_BAR: Display a 'Progress Bar'
% 
% Useage: progress_bar('init',progress,title)
%         Initialises the progress bar.  Input 'progress' is
%         optional, assumed zero, otherwise a decimal
%         value between 0:1.  Input 'title' is optional,
%         but can be used to indicate the process of the
%         progress bar.
%
% Useage: progress_bar('Set',progress)
%         Updates the progress bar.  Input 'progress' is
%         a decimal value between 0:1.
%
% Useage: progress_bar('Clear')
%         Clears the progress bar.
%

%-----------------------------------------------------------------------
% @(#)spm_progress_bar.m	2.1 John Ashburner 99/05/17
% Modified 03/2002, Darren.Weber@flinders.edu.au
%                   - removed spm specific references
%                   - modified inputs/input handling to my liking
%                   - progress bar is now a patch object and the
%                     erasemode is xor with no figure backingstore
% Modified 03/14/04, Vasco Curdia

%-----------------------------------------------------------------------
if isequal(nargin,0),
	action = 'init';
end

action = lower(action);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialize
if strcmp(action,'init'),
    global tim
    tim = clock;
    timestr = sprintf('  Began %02.0f:%02.0f:%02.0f',tim(4),tim(5),tim(6));
	
    if exist('title','var'),
        if ~isempty(title),
            name = strcat(sprintf('%s  -',title),timestr);
        else
            name = strcat('Progress  -',timestr);
        end
    else
        name = strcat('Progress  -',timestr);
    end
    
    fg = figure('MenuBar','none',...
                'NumberTitle','off',...
                'Name',name,...
                'Tag','ProgressBar',...
                'BackingStore','off',...
                'Pointer','watch',...
                'Position',[1 1 500 60]);
    movegui(fg,'northeast');
    ax = axes('Position',[0.1 0.35 0.8 0.2],...
              'YTick',[],...
              'Xlim', [0 1],'Ylim',[0 1],...
              'Box','on',...
              'FontSize',8,...
              'Tag','ProgressBarAxis',...
              'Parent',fg);
    
    xlab = get(ax,'xticklabel');
    xlab = str2num(xlab) * 100;
    xlab = num2str(xlab);
    set(ax,'xticklabel',xlab)
    
    if exist('progress','var'),
        if ~isempty(progress), setpb(fg,progress);
        else
            progress = 0;
            setpb(fg,progress);
        end
    else
        progress = 0;
        setpb(fg,progress);
    end
    
    drawnow;
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set
elseif strcmp(action,'set'),
    
    fg = findobj('Tag','ProgressBar');
    
    if ~isempty(fg), setpb(fg,progress); end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% clear
elseif strcmp(action,'clear'),
    fg = findobj('Tag','ProgressBar');
    if ~isempty(fg), close(fg);	end;
end;

return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function setpb(fig,progress)
    
    pbaxis = findobj(fig,'Tag','ProgressBarAxis');
    
    if ~isempty(pbaxis),
        
        vert = [0 0; progress 0; progress 1; 0 1];
        face = [1 2 3 4];
        
        pbpatch = findobj(fig,'Tag','ProgressBarPatch');
        
        if ~isempty(pbpatch),
            set(pbpatch,'Vertices',vert);
        else
            pbpatch = patch('Faces',face,'Vertices',vert,'FaceColor',[54/255 95/255 114/255],...
                'Tag','ProgressBarPatch',...
                'EraseMode','none',...
                'Parent',pbaxis);
        end
        
        global tim
        vcet=etime(clock,tim)/(60^2);
        vch=fix(vcet);
        vcet=(vcet-vch)*60;
        vcm=fix(vcet);
        vcet=(vcet-vcm)*60;
        vcs=fix(vcet);
        vcr=(vcet-fix(vcs))*100;
        vcet=etime(clock,tim)*(1-progress)/progress/(60^2);
        vchr=fix(vcet);
        vcet=(vcet-vchr)*60;
        vcmr=fix(vcet);
        vcet=(vcet-vcmr)*60;
        vcsr=fix(vcet);
        vcrr=(vcet-fix(vcsr))*100;
        title = get(pbaxis,'Title');
        set(title,'string',sprintf('%5.0f%% Complete  Time Elapsed: %.0fh%.0fm%.0fs%.0f  Time Remaining: %.0fh%.0fm%.0fs%.0f',...
                fix(100*progress),vch,vcm,vcs,vcr,vchr,vcmr,vcsr,vcrr),'EraseMode','xor');
        drawnow;
        
        figure(fig);
    end;
return

