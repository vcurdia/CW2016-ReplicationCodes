function vcpath(RootFolder)

% vcpath
%
% Add MatlabCodes folders to path
%
% ...........................................................................
%
% Created: January 26, 2016 by Vasco Curdia
% 
% Copyright (C) 2016-2017 Vasco Curdia

%% -------------------------------------------------------------------

if ~exist('RootFolder','var'), RootFolder = strrep(pwd,'VC_Tools',''); end

MyPathDependencies = {...
    'VC_DSGE';
    'VC_BayesianEstimation';
%     'VC_CREstimation';
    'VC_LQ';
    'VC_Tools';
%     'RedsSolds';
    'Sims_gensys';
    'Sims_KF';
    'Sims_optimize';
    'Sims_VARtools';
    'JaeWon';
    'Misc';
    };

fprintf('Adding folders to path:\n');
for j=1:length(MyPathDependencies)
    pj = [RootFolder,MyPathDependencies{j}];
    disp(pj)
    addpath(pj)
%     MyPathDependencies{j} = [RootFolder,MyPathDependencies{j}];
end
% addpath(MyPathDependencies{:})

%% -------------------------------------------------------------------
