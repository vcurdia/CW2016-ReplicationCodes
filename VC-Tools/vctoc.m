function out=vctoc(t)

% vctoc
%
% This file analyzes and reports elapsed time in more detailed fashion, as
% compared to simple toc command.
% NOTE: it assumes that toc command was used before.
%
% Usage:
%
%   vctoc
%   uses global time variable toc as reference
%
%   vctoc(t)
%   uses toc(t), where tic was assigned to t previously or t is already the
%   time elapsed
%
%   out = vctoc
%   out = vctoc(t)
%   same as above but saves txt in out, instead of displaying it
%
% .........................................................................
%
% Created: October 19, 2002 by Vasco Curdia
% 
% Copyright 2002-2017 by Vasco Curdia

% -------------------------------------------------------------------------

if nargin==1
    if strcmp(class(t),'double')
        sElapsed = t;
    else
        sElapsed = toc(t);
    end
else
    sElapsed = toc;
end

hElapsed = fix(sElapsed/(60^2));
sElapsed = sElapsed - hElapsed*(60^2);
mElapsed = fix(sElapsed/60);
sElapsed = sElapsed - mElapsed*60;
rElapsed = sElapsed - fix(sElapsed);
sElapsed = sElapsed - rElapsed;
rElapsed = fix(100*rElapsed);

txt = 'Elapsed time is';
if hElapsed>0, txt = sprintf('%s %.0fh',txt,hElapsed); end
if (mElapsed>0)||(hElapsed>0), txt = sprintf('%s %.0fm',txt,mElapsed); end
txt = sprintf('%s %.0fs %02.0f',txt,sElapsed,rElapsed);

if nargout==1
    out = txt;
else
    disp(txt)
end
