function y=iif(varargin)
% iif
% 
% Inline Conditional. Inspired by Mathworks blog post on Introduction to 
% Functional Programming with Anonymous Functions, Part I
% http://blogs.mathworks.com/loren/2013/01/10/introduction-to-functional-programming-with-anonymous-functions-part-1/
%
% Unlike that example, this version does not allow multiple output arguments, 
% and simply evaluates expressions in caller workspace.
%
% Usage:
% y = iif( if this,      then run this, ...
%          else if this, then run this, ...
%          ...
%          else,         then run this );
%
% Example:
%   gam = 2;
%   eta = 0.5;
%   z = iif(eta==1, 'gam', true, 'gam^(1/(1-eta))' )
%   z =
%        4
%   eta = 1;
%   z = iif(eta==1, 'gam', true, 'gam^(1/(1-eta))' )
%   z =
%        2
%
% ..............................................................................
%
% Created: October 5, 2016 by Vasco Curdia
% 
% Copyright 2016 by Vasco Curdia

% Original  anonymous function code:
% iif = @(varargin) varargin{2*find([varargin{1:2:end}], 1, 'first')}();

% The following works but may be very inefficient if iif called often in same 
% workspace because each time it is creating a new liif function
% liif = @(varargin) varargin{2*find([varargin{1:2:end}], 1, 'first')}();
% varargout = cell(1,nargout);
% [varargout{:}] = liif(varargin{:});

% this is more basic but efficient for simple calls
y = evalin('caller',varargin{2*find([varargin{1:2:end}], 1, 'first')});

