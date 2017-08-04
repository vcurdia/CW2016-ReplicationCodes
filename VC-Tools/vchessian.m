function [H,badg] = vchessian(fcn,x,varargin)

% vchessian
%
% Computes a numerical Hessian for function "fcn"
%
%   [H,badg] = vchessian(fcn,x,varargin)
%
% inputs:
%
%	fcn
%   the name of the function to be differentiated (can be a vector of
%   functions)
%
%   x
%   vector of inputs for function around which differentiated is performed
%
%   varargin (optional)
%   optional additional parameters to be used in the function
%
% outputs:
%   H
%   Hessian matrix of the function
%   
%   badg (optional)
%   Flag for bad gradient: 0 if normal; 1 if bad gradient found.
%
% Required mfiles
%   vcnumjacobian
%
% ..............................................................................
%
% Created: August 24, 2008 by Vasco Curdia
% 
% Copyright 2008-2017 by Vasco Curdia

%% ------------------------------------------------------------------------

[H,badg] = vcnumjacobian(@(xx)vcnumjacobian(fcn,xx,varargin{:}),x);

%% ------------------------------------------------------------------------

