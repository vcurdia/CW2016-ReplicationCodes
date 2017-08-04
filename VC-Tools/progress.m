StartTime = clock;
h = waitbar(0,'Please wait...');
N = 30;
%rect = [left, bottom, width, height]
rect = [470.2500  504.0000  270.0000   56.2500];
set(h,'Position',rect);

for i=1:N, % computation here %
    timeelapsed = num2str(etime(clock,StartTime),'%6.2f');
    timeremain = num2str(etime(clock,StartTime)/(i/N)-etime(clock,StartTime),'%6.2f');
    text = ['\bf' num2str(i/N*100,'%6.2f') '\rm% complete: \bf' timeelapsed '\rm sec elapsed / \bf' timeremain '\rm sec remain'];
    waitbar(i/N,h,text)
    pause(1);
    
end
close(h) 


%% -------------------
disp (sprintf('Run Time =            %6.2f sec', etime(clock,StartTime)));
disp (datestr(now))