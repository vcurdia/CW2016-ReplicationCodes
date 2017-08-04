function [G1,C,impact,fmat,fwt,ywt,gev,eu]=vcgensys(...
  Gamma0,Gamma1,Const,Psi,Pi,Author,varargin)

% vcgensys
%
% Allows use of alternative versions of gensys
%
% Option: Author (string)
% If set to 'CS' (default) it uses the original verion. If set to 'JW' it uses
% the fast gensys from Jae Won and if it does not yield a normal solution it
% runs the original gensys.
%
% ..............................................................................
% 
% Created: February 14, 2011 by Vasco Curdia
% Updated: July 26, 2011 by Vasco Curdia
%          Allow for the ommission of div, because the default value was causing
%          problems.
% 
% Copyright 2011 by Vasco Curdia

%% -----------------------------------------------------------------------------

%% Preamble
eu=[0;0];
if ~exist('Author','var'), Author = 'CS'; end

%% Run JW
if strcmp(Author,'JW')
  [G1,C,impact,fmat,fwt,ywt,gev,eu] = ...
    fastgensysJaeWonvb(Gamma0,Gamma1,Const,Psi,Pi,varargin{:});
end

%% Run CS
if all(eu(:)==1), return, end
[G1,C,impact,fmat,fwt,ywt,gev,eu] = ...
  gensysvb(Gamma0,Gamma1,Const,Psi,Pi,varargin{:});

%% -----------------------------------------------------------------------------