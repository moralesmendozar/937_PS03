function [nextst,outstr]=calcStateTransition(obj,point,solvec,mode,varargin)
params=obj.Params;
% unpack relevant parameters
alpha=params.alpha;
Kbar=params.Kbar;
phi=params.phi;

% extract state variables
exst=point(1);
kF=point(2);
levF=point(3);
a=obj.Exogenv.pts_all(exst,1);
alower=obj.Exogenv.pts_all(exst,2);

% compute next period's states
if isempty(varargin)
    State_next=obj.Tfct.evaluateAt(point);
else
    State_next=varargin{1};
end

% extract endogeous variables
%R=exp(solvec(1));
q=exp(solvec(2));
cF=exp(solvec(3));
cG=exp(solvec(4));
kGpol=exp(solvec(5));
%muF=solvec(6);

% market clearing for capital
kG = Kbar - kF;
kFpol = Kbar - kGpol;
% aggregate output
Y = a*kF + alower*kG^alpha;
% farmer wealth
RbF = levF * q*kF;
wF = (a+q)*kF - RbF;
% wealth distribution
wG = alower*kG^alpha + q*kG + RbF;

% bond portfolio choices
bFpol = wF - cF - q*kFpol - phi/2*(kFpol/kF - 1)^2 * kF;
bGpol = wG - cG - q*kGpol;

if mode>0
    % simulation, mode contains number of next period's state
    exst=mode;
    exnpt=size(obj.Exogenv.pts_all,1);
    cind=obj.Exogenv.pts_all(exst,end);
    wFnext=State_next(exst,1);
    kFnext=State_next(exnpt+exst,1);
    nextst=[cind,wFnext,kFnext];
else
    % solution mode, compute next period's state variables for all
    % possible Markov states
    cind=obj.Exogenv.pts_all(:,end);
    nst=obj.Exogenv.exnpt;
    wFnext=State_next(1:nst,1);
    kFnext=State_next(nst+(1:nst),1);
    % matrix of next period states
    nextst=[cind,wFnext,kFnext];
end

addvars=struct('bFpol',bFpol,...
    'bGpol',bGpol,...
    'kGpol',kGpol,...
    'kF',kF,...
	'levF',levF,...
    'Y',Y,...
    'wG',wG,...
    'wF',wF);

outstr=struct;
outstr.addvars=addvars;
outstr.exstvec=[a;alower];

end