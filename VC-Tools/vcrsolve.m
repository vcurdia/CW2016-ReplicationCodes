function varargout = vcrsolve(sys)

% this function solves a system by recursive substitution
% NOTE: it allows only for unique solutions, otherwise it searches for
%       complex solutions and eliminates them. If there are multiple real
%       solutions then it yields an error message.
%
% this function requires the symbolic toolbox
%
% ..............................................................................
%
% Created: September 3, 2004 by Vasco Curdia
% Updated: July 26, 2011 by Vasco Curdia
% 
% Copyright 2004-2011 by Vasco Curdia


%--------------------------------------------------------------------------

% find variables and create a vector of them
vars = findsym(sys);
vars(findstr(',',vars))='';
eval(['syms ',vars])
vars = eval(['[',vars,']']);
vars_sol = eval(vars);
nvars = length(vars);
for j=1:nvars
    var_names{j}=char(vars(j));
end

% remaining number of variables to be solved for
rvars = vars;
nrvars = length(rvars)
rsys = sys;
nrsys = length(rsys);

% initialize error checks
constant = 0;

% keep recursive iterations until all variables are solved for
while nrvars>0
    if nrvars>nrsys
        error('more unknowns than equations')
    end
    solved = 0;
    for j=1:nrsys
        if solved, break, end
        veq = eval(['[',findsym(rsys(j)),']']);
        if length(veq)==1
            % if equation contains only one variable then solve for it
            % test whether it is not a valid equation 
            if simple(rsys(j))~=0
                if isempty(findstr('+',char(simple(rsys(j)))))
                    if isempty(findstr('+',char(-simple(rsys(j)))))
                        eval([char(veq),'=sym(0);'])
                    else
                        eval([char(veq),'=solve(rsys(j));'])
                    end
                else
                    eval([char(veq),'=solve(rsys(j));'])
                end
                nsol = length(eval(veq));
                for js=nsol:-1:2
                    vn=eval(veq);
                    if vn(js)==vn(1)
                        eval([char(veq),'(js)=[];'])
                    elseif round(1e+12*double(subs(rsys(j),char(veq),...
                            vn(js))))*1e-12~=0
                        % if for some reason this solution is not
                        % correct then ignore it
                        eval([char(veq),'(js)=[];'])
                    elseif ~isreal(vn(js))
                        % if this is not a real solution, ignore it
                        eval([char(veq),'(js)=[];'])
                    elseif isreal(vn(js))&~isreal(vn(1))
                        % if it is real and the first one is not,
                        % replace it
                        vn(1)=vn(js);
                        eval([char(veq),'(1)=',char(veq),'(js);'])
                        eval([char(veq),'(js)=[];'])
                    else
                        fprintf('equation %.0f',j)
                        rsys(j)
                        eval(char(veq))
                        error(['multiple solutions found for at least',...
                                ' one variable'])
                    end
                end
                eval(char(veq))
            else
                rsys(j)=[];
            end
            % leads the system to simplify before moving to other 
            % equations
            solved = 1;
        end
    end
    
    % if no equations were found with only one variable then proceed
    % iteratively by searching for 2,3,...,nrvars variables and start
    % solving for one and plug it in and so on
    for nv=2:nrvars
        for j=1:nrsys
            if solved, break, end
            veq = eval(['[',findsym(rsys(j)),']']);
            if length(veq)==nv
                % test whether it is a valid equation
                if simple(rsys(j))==0
                    rsys(j)=[];
                    solved=1;
                end
                
                % if equation contains only nv variables then solve for
                % the first one as depending on the others unless it yields
                % multiple solutions
                for jj=1:nv
                    if solved, break, end
                    soltmp=solve(rsys(j),veq(jj));
                    if length(soltmp)==1
                        eval([char(veq(jj)),'=soltmp'])
                        solved=1;
                        % then lead the system to simplify before moving to
                        % further stages of search for variables
                        break
                    end
                    clear soltmp
                end
                
                % if after searching for a single solution in all variables
                % present in this equation no such thing happened, move to
                % the next equation and if not any for this, go to more
                % variables per equation and keep going...
            end
        end
    end
    
    if ~solved
        % if at this stage there was no progress in simplifying the model
        % then this means that multiple solutions were found at an
        % intermediate step
        nrsys
        nrvars
        rsys
        error(['multiple solutions found in the process of simplifying',...
                'the system.'])
    end
    
    vars = eval(vars);
    rsys = eval(rsys);
    idxsys = [];
    for jj=1:length(rsys)
        if isempty(findsym(rsys(jj)))
            if round(double(rsys(jj))/1e-12)*1e-12==0
                idxsys=[idxsys,jj];
            else
                constant = 1;
                idxsys=[idxsys,jj];
            end
        elseif isempty(findstr(',',findsym(rsys(jj))))
            % in the case of only one symbolic variable check whether it
            % simplifies to 0=0 and if so eliminate this equation
            if simple(rsys(jj))==0
                idxsys=[idxsys,jj];
            else
                if isempty(findstr('+',char(rsys(jj))))
                    if isempty(findstr('+',char(-rsys(jj))))
                        eqres=subs(rsys(jj),findsym(rsys(jj)),1);
                        if round(eqres/1e-12)*1e-12==0
                            idxsys=[idxsys,jj];
                        end
                    end
                end
            end
        end
    end
    rsys(idxsys)=[]
    if ~isempty(rsys)
        rvars = eval(['[',findsym(rsys),']']);
        nrvars = length(rvars)
        nrsys = length(rsys);
    else
        nrvars = 0;
    end
end

% evaluate variables
for jj = nvars
    for j=1:nvars
        eval([var_names{j},'=',char(vars(j)),';'])
    end
    vars = eval(vars);
end

% show warning message if error occurred
if constant
    warning(['Warning - Possible wrong solution. Solution does not',...
        ' satisfy equations:'])
    ssys = double(eval(sys))
    idx = find(round(double(ssys)/1e-12)*1e-12~=0)
    sys(idx)
end

% prepare the output as a structure or as a cell array
if nargout<=1
    for j=1:nvars
        eval(['solvesol.',var_names{j},' = vars(j);'])
	end
    varargout{1} = solvesol;
elseif nargout==nvars
	for j=1:nvars
        varargout{j} = vars(j);
	end
else
    error([int2str(nvars) ' variables does not match ' ...
             int2str(nargout) ' outputs.'])
end

