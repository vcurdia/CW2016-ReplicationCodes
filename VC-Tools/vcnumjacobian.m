function [g,badg] = vcnumjacobian(fcn,x,varargin)

% Computes a numerical jacobian for function "fcn"
%
%   [g,badg] = vcnumjacobian(fcn,x,varargin)
%
% NOTE: this is a generalization of Chris Sims' code numgrad.m
%
% Inputs:
%  fcn
%  the name of the function to be differentiated (can be a vector of 
%  functions)
%  
%  x
%  vector of inputs for function around which differentiated is performed
%
%  varargin
%  optional additional parameters to be used in the function
%
% Outputs:
%
%  g
%  The vector/matrix with the facobian of the function.
%  If it is a single function then it is a (nx x 1) matrix.
%  If it is a vector of fucntions then it is a (nf x nx) matrix
%
%  badg (optional)
%  Flag for bad gradient: 0 if normal; 1 if bad gradient found.
%
% .........................................................................
%
% Created: June 16, 2006 by Vasco Curdia
% 
% Copyright 2006-2017 by Vasco Curdia

%% ------------------------------------------------------------------------

delta = 1e-6; % original 1e-6
nx = length(x);
tvec=delta*eye(nx);
f0 = feval(fcn,x,varargin{:});
nf = length(f0);
g=zeros(nf,nx);
badg=0;
for jx=1:nx
    if size(x,1)>size(x,2)
        tvecv=tvec(jx,:);
    else
        tvecv=tvec(:,jx);
    end 
    g0 = (feval(fcn,x+tvecv',varargin{:})-f0)/delta;
    if abs(g0)< 1e15
        g(:,jx)=g0; % good gradient
    else
        warning('bad gradient')
        g(find(~(abs(g0)< 1e15)),jx)=0;
        badg=1;
        % return
        % can return here to save time if the gradient will never be
        % used when badg returns as true.
    end
end
if nf == 1, g = g'; end
