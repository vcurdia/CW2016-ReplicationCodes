% PlotAll
%
% Commands the execution of all plots of all models under all 
% specifications and policy rules, including comparison of different models
% and policies
%
% .........................................................................
%
% Copyright 2008-2009 by Vasco Curdia and Michael Woodford
% Created: February 26, 2008
% Updated: February 26, 2009

%% ------------------------------------------------------------------------

%% Preamble

clear all
tic
setpath


if ~isempty(dir('IRF*.eps'))
    mkdir('IRFTEMP')
    movefile('IRF*.eps','IRFTEMP')
    fprintf('\nWARNING: Found files IRF*.eps in current directory.')
    fprintf('\n         They were moved to directory: IRFTEMP\n\n')
end

Models = {'FF','NoFF','RepHH'};

% nCP = 1; Policy{nCP} = {'LQ','PiStab','Taylor'};
nCP = 1; Policy{nCP} = {'LQ','PiStab','Taylor','FlexTarget','TaylorYn'};

Pers = 1;
End = 0:1;
NoRes = 0;
NoSpread = 0;
NoDist = 0;
GDebt = 0;
SmSigma = 0:1;
LowSigma = 0;

nM = length(Models);
nP = length(Policy{1});

%% set the specifications to consider
jE = 0;
for isPers=Pers
    for isEnd=End
        for isNoRes=NoRes
            for isNoSpread=NoSpread
                for isNoDist=NoDist
                    for isGDebt=GDebt
                        for isSmSigma=SmSigma
                            for isLowSigma=LowSigma
                                if isPers
                                    FileNameSuffix = '_Pers';
                                else
                                    FileNameSuffix = '_NoPers';
                                end
                                if isEnd
                                    FileNameSuffix = sprintf('%s_End',FileNameSuffix);
                                else
                                    FileNameSuffix = sprintf('%s_Exo',FileNameSuffix);
                                end
                                if isNoRes, FileNameSuffix = sprintf('%s_NoRes',FileNameSuffix); end
                                if isNoSpread, FileNameSuffix = sprintf('%s_NoSpread',FileNameSuffix); end
                                if isNoDist, FileNameSuffix = sprintf('%s_NoDist',FileNameSuffix); end
                                if isGDebt, FileNameSuffix = sprintf('%s_GDebt',FileNameSuffix); end
                                if isSmSigma, FileNameSuffix = sprintf('%s_SmSigma',FileNameSuffix); end
                                if isLowSigma, FileNameSuffix = sprintf('%s_LowSigma',FileNameSuffix); end
                                jE = jE+1;
                                ExerciseName{jE} = FileNameSuffix;
                            end
                        end
                    end
                end
            end
        end
    end
end
nE = jE;

%% ------------------------------------------------------------------------

%% For each exercise...
for jE=1:nE
    fprintf('\n****%s****',repmat('*',1,length(ExerciseName{jE})))
    fprintf('\n*   %s   *',ExerciseName{jE})
    fprintf('\n****%s****\n\n',repmat('*',1,length(ExerciseName{jE})))

%% Single plot for all models and policies
% fprintf('\nSingle plot for each model and policy:')
% fprintf('\n======================================\n')
% for jM=1:1%nM %do it only for FF, we don't really care about the others...
%     fprintf('\nModel: %s\n',Models{jM})
%     for jP=1:nP
%         fprintf('\nPolicy: %s\n',Policy{1}{jP})
%         IRFPlotCompare(Models(jM),Policy{1}(jP),ExerciseName{jE})
%         close all
%     end
% end

%% Multiple models, single policy
fprintf('\nMultiple models, single policy...\n')
% fprintf('\n===============================\n')
for jP=1:nP
%     fprintf('\nPolicy: %s\n',Policy{1}{jP})
    IRFPlotCompare(Models,Policy{1}(jP),ExerciseName{jE})
    close all
end

%% Single model, multiple policies
fprintf('Single model, multiple policies...\n')
% fprintf('\n================================\n')
for jM=1%:nM
%     fprintf('\nModel: %s\n',Models{jM})
    for jCP=1:nCP
        if ~(jCP>1&&~strcmp(Models(jM),'FF'))
            if any(ismember({'TaylorSP50Bm50','TaylorBSPm50'},Policy{jCP}))...
                    && ~isempty(strfind(ExerciseName{jE},'SmSigma'))
                continue
            end
            IRFPlotCompare(Models(jM),Policy{jCP},ExerciseName{jE})
            close all
        end
    end
end

%% move files and compile
fprintf('\nMoving and compiling files...\n\n')
dirName = strrep(ExerciseName{jE},'_','');
mkdir(dirName)
movefile('IRF*.eps',dirName)
% eval(['! c:\eps2pdf\eps2pdf /d=',pwd,'\',dirName])
cd(dirName)
epstopdfall
cd ..
delete([dirName,'\*.eps'])

%% end exercises
end

%% ------------------------------------------------------------------------

disp(' '),vctoc,disp(' ')
