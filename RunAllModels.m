% RunAllModels
%
% Commands the execution of all models under all specifications
%
% .........................................................................
%
% Copyright 2008-2009 by Vasco Curdia and Michael Woodford
% Created: February 26, 2008
% Updated: February 27, 2009

%% ------------------------------------------------------------------------

%% Run all experiments
clear all
tic
setpath


%% Settings
Models = {'FF','NoFF','RepHH'};

isPers = 1;
isNatVars = 1;

End = 0:1;
NoRes = 0;
NoSpread = 0;
NoDist = 0;
GDebt = 0;
SmSigma = 0:1;
LowSigma = 0;

nM = length(Models);
nExercise = 0;
for jM=1:nM
    for isEnd=End
        for isNoRes=NoRes
            for isNoSpread=NoSpread
                for isNoDist=NoDist
                    for isGDebt=GDebt
                        for isSmSigma=SmSigma
                            for isLowSigma=LowSigma
                                nExercise = nExercise+1;
                                Options = {};
                                if ~isPers, Options{end+1} = 'NoPers'; end
                                if isNatVars, Options{end+1} = 'NatVars'; end
                                if isEnd
                                    Options{end+1} = 'End';
                                else
                                    Options{end+1} = 'Exo';
                                end
                                if isNoRes, Options{end+1} = 'NoRes'; end
                                if isNoSpread, Options{end+1} = 'NoSpread'; end
                                if isNoDist, Options{end+1} = 'NoDist'; end
                                if isGDebt, Options{end+1} = 'GDebt'; end
                                if isSmSigma, Options{end+1} = 'SmSigma'; end
                                if isLowSigma, Options{end+1} = 'LowSigma'; end
                                ExerciseCmd{nExercise} = sprintf('IntModel%s',Models{jM});
                                ExerciseOptions{nExercise} = Options;
                            end
                        end
                    end
                end
            end
        end
    end
end
for jE=1:nExercise
	feval(ExerciseCmd{jE},ExerciseOptions{jE}{:})
end

%% ------------------------------------------------------------------------

%% Elapsed time
fprintf('\n\n%s\n\n',vctoc)

%% ------------------------------------------------------------------------
