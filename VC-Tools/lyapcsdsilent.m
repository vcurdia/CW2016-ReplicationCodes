function  [F,eno]=lyapcsdsilent(R,S,Silent)
%function  [F,eno]=lyapcsd(R,S)
% solves RFR'+S=F
% Will fail if R has a root that is one in absolute value.
if nargin<3,Silent=0;end
IR=eye(size(R))+R;
[F,eno]=lyapcssilent(-inv(IR),R'/IR',(IR\S)/IR',Silent);
