function [btTn,StTn] = distsmoother(btt,Stt,bn,A,omega)

% distsmoother
%
% Disturbance smoother recursion.
%
% Usage:
%   [btT,StT]=distsmoother(btt,Stt,bn,A,omega)
%
% State evolution equation is
%   bn=A*bt+e,  Var(e)=omega
%
% with
%   bt|t ~ N(btt,Stt) -- from Kalman Filter
%   bt|T,bn ~ N(btTn,StTn)
% btt and bn are assumed to be loaded as row vectors and btTn is also in
% the same format
%
% ..............................................................................
%
% Created: October 29, 2007 by Vasco Curdia
% 
% Copyright 2007-2017 by Vasco Curdia

%% ------------------------------------------------------------------------

AS = A*Stt;
G = AS*A'+omega;

%****************************
% SAGI=AS'/G;
% line below may slightly slow the routine, but makes it robust vs.
% cases where part or all of the state vector is known with certainty.
% SAGI = AS'*pinv(G);
%*****************************

btTn = btt+(AS'*(G\(bn'-A*btt')))';
StTn = Stt-AS'*(G\AS);

%% ------------------------------------------------------------------------
