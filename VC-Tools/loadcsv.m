function Data = loadcsv(DataFile)

% loadcsv
%
% Loads a csv file and separates dates, variables and data matrix.
%
% Syntax:
%
%   Data = loadcsv(DataFile)
%
% .........................................................................
%
% Created: May 18, 2015 by Vasco Curdia
%
% Copyright 2015-2017 by Vasco Curdia

%% ------------------------------------------------------------------------

Data = importdata(DataFile);
Data.TimeIdx = Data.textdata(2:end,1)';
Data.TimeStart = Data.TimeIdx{1};
Data.TimeEnd = Data.TimeIdx{end};
Data.T = length(Data.TimeIdx);
Data.Var = {Data.textdata{1,2:end}};
Data.NVar = length(Data.Var);
Data.Values = [Data.data;nan(Data.T-size(Data.data,1),Data.NVar)];
Data = rmfield(Data,{'textdata','data'});


%% ------------------------------------------------------------------------
