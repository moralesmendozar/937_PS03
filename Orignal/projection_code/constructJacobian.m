function Jac = constructJacobian(obj)

params = obj.Params;
% unpack relevant parameters
alpha=params.alpha;
Kbar=params.Kbar;
phi=params.phi;
betaF=params.betaF;
betaG=params.betaG;
alpha=params.alpha;
Kbar=params.Kbar;
gamma=params.gamma;
theta=params.theta;
phi=params.phi;

% Define symbolic state variables
syms a alower kF levF;

% Define symbolic policy functions
syms R q cF cG kGpol muFplus muFminus lamFplus lamFminus lamGplus lamGminus;

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

dim = [1,obj.Exogenv.exnpt];
cF_next = sym('cF_next',dim(end:-1:1));
cG_next = sym('cG_next',dim(end:-1:1));
q_next = sym('q_next',dim(end:-1:1));
kFpol_next = sym('kFpol_next',dim(end:-1:1));
prnext = sym('prnext',dim);

if params.use_min_q
	syms min_q
else
	min_q = q;
end

a_next=obj.Exogenv.pts_all(:,1);
alower_next=obj.Exogenv.pts_all(:,2);

% Marginal values for gatherer
MUg_next = betaG * (cG_next).^(-gamma);
exp_EEG_b = prnext*MUg_next;
exp_EEG_k = prnext*(MUg_next.*(alower_next*alpha.*kGpol.^(alpha-1) + q_next));

% Marginal values for farmer
MUf_next = betaF * (cF_next).^(-gamma);
exp_EEF_b =  prnext*MUf_next;
exp_EEF_k= prnext*(MUf_next.*(a_next + q_next + phi/2*((kFpol_next./kFpol).^2 - 1)));

% FOCs for gatherer
MUg = cG^(-gamma);
fx(1) = 1/R - exp_EEG_b / MUg;
fx(2) = q - lamGplus - exp_EEG_k / MUg;
fx(3) = kGpol - lamGminus;


% FOCs for farmer
MUf = cF^(-gamma);
fx(4)= 1 - muFplus - R * exp_EEF_b / MUf;
fx(5)= q + phi*(kFpol/kF - 1) - lamFplus - theta*min_q*muFplus - exp_EEF_k / MUf;
fx(6)= theta*min_q*kFpol + bFpol - muFminus;
fx(7)= kFpol - lamFminus;


% market clearing for bonds
fx(8)= bGpol + bFpol;

symJac = jacobian(fx,[R q cF cG kGpol muFplus muFminus lamFplus lamFminus lamGplus lamGminus]);
symJac = simplify(symJac);

args = [R q cF cG kGpol muFplus muFminus lamFplus lamFminus lamGplus lamGminus ...
	a alower kF levF ...
	reshape(cF_next,dim) reshape(cG_next,dim) reshape(q_next,dim) ...
	reshape(kFpol_next,dim) prnext ];

if params.use_min_q
	args = [args, min_q];
end

Jac = matlabFunction(symJac,'Vars', args, ...
	'File','jacfun.m');

end