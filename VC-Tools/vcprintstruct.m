function []=vcprintstruct(StructName,addText)

% Prints each variable of the structure according to its name in the
% structure 
%
% []=vcprintstruct(StructName,addText)
%
% inputs:
%   StructName - a string with the name of a structure with at least one 
%                field called "name", containing the name of the variables 
%                to show (and these variables must corresponde to variables 
%                defined in the base environment as double).
%   addText -   (optional) adds this text to the variable names.
%
% ..............................................................................
%
% Created: August 26, 2005 by Vasco Curdia
% Updated: July 26, 2011 by Vasco Curdia
% 
% Copyright 2005-2011 by Vasco Curdia

%% ------------------------------------------------------------------------

if nargin == 1
    addText = [];
end
vStruct = evalin('caller',StructName);
nV = length(vStruct);
if nV == 1
    disp(sprintf([StructName,' = {''',vStruct(1).name,''', %f};'],double(evalin('caller',[vStruct(1).name,addText]))))
else
    nGap = length(StructName)+4;
    disp(sprintf([StructName,' = {''',vStruct(1).name,''', %f;'],double(evalin('caller',[vStruct(1).name,addText])))) 
    for j=2:nV-1
        disp(sprintf([repmat(' ',1,nGap),'''',vStruct(j).name,''', %f;'],double(evalin('caller',[vStruct(j).name,addText]))))
    end
    disp(sprintf([repmat(' ',1,nGap),'''',vStruct(nV).name,''', %f};'],double(evalin('caller',[vStruct(nV).name,addText]))))
end
