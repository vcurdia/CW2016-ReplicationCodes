function ys=seasadj(y,freq,tsize)

% Conversion of GAUSS code for applying the X11 filter
%
% GAUSS original from M. Watson
% seasadj.prc, 4/11/00, mww
%   Seasonal Adjustment
%   Seasonally adjust a series using a X11 if the series
%   fails a "non-seasonal" pretest
%   Adjustment fails if
%   (i) too few obs for pretest
%   (ii) missing values in the "middle" of the series
%   When adjustment fails then scalar missing value is returned
%   Inputs:
%   Y -- Series to be adjusted
%   freq -- m - monthly; q - quarterly
%   tsize -- size of pretest
%   The returns are:
%   YS == seasonally adjusted Y (scalar missing value if proc fails)
%
% ..............................................................................
%
% Created: October 29, 2002 by Vasco Curdia
% Updated: July 26, 2011 by Vasco Curdia
% 
% Copyright 2002-2011 by Vasco Curdia

ys=y;

if ~isempty(tsize)
    san=sadet(ys,freq,tsize); %checks for seasonality and applies x11 as required
    if san==0
        disp('No adjustment necessary based on pretest. Ysa=Y');
        return
    elseif san==2
        error('Too few obs for pretest');
    end
end

temp=[ys [1:1:size(ys,1)]'];
temp(any(isnan(temp)'),:)=[];
ys=temp(:,1);
t1=temp(1,2);           % Index of First Non-missing value
t2=temp(size(ys,1),2);    % Index of Last Non-missing value
t2a=t1-1+size(ys,1);      % Val of t2 if no internal Missing
if t2~=t2a
    error('Missig values in interior of series');
end
ys=x11ar(ys,freq);
if t1~=1
    ys=[repmat(NaN,t1-1,1);ys];
end
if t2~=size(y,1)
    ys=[ys;repmat(NaN,size(y,1)-t2,1)];
end
    
