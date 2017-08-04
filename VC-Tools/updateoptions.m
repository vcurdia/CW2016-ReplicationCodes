function op = updateoptions(op,varargin)

% updateoptions
%
% updates a structure of options using either a partially filled structure of 
% options, or a cell array, or both
%
% ...........................................................................
%
% Created: December 20, 2016 by Vasco Curdia
% 
% Copyright (C) 2016-2017 Vasco Curdia

opfields = fieldnames(op);

if length(varargin)>0 && isstruct(varargin{1})
    opnew = varargin{1};
    opnewfields = fieldnames(opnew);
    for jop=1:length(opnewfields)
        opj = opnewfields{jop};
        if ismember(opj,opfields)
            op.(opj) = opjupdate(op.(opj),opnew.(opj));
        else
            op.(opj) = opnew.(opj);
        end
%         if ismember(opj,opfields)
%             op.(opj) = opnew.(opj);
%         else
%             error('Option name not recognized: %s',opj)
%         end
    end
    varargin(1) = [];
end

for jop=1:(length(varargin)/2)
    opj = varargin{(jop-1)*2+1};
    if ismember(opj,opfields)
        op.(opj) = opjupdate(op.(opj),varargin{jop*2});
    else
        op.(opj) = varargin{jop*2};
    end
%     if ismember(opj,opfields)
%         op.(opj) = varargin{jop*2};
%     else
%         error('Option name not recognized: %s',opj)
%     end
end

end

function opnew = opjupdate(op,opnew)
    if isstruct(op)
        opnew = updateoptions(op,opnew);
    end
end
