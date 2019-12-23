function [fx,J,V]=calcEquations(obj,exst,nextst,solvec,instr,mode,varargin)

% allocate result
fx=zeros(obj.NSOL,1);
params=obj.Params;
% unpack relevant parameters
alpha=params.alpha;
Kbar=params.Kbar;
gamma=params.gamma;
theta=params.theta;
phi=params.phi;
betaF=params.betaF;
betaG=params.betaG;

% extract endogeous variables
R=exp(solvec(1));
q=exp(solvec(2));
cF=exp(solvec(3));
cG=exp(solvec(4));
kGpol=exp(solvec(5));
muF=solvec(6);
lamF=solvec(7);
lamG=solvec(8);

% multiplier transformations
muFplus=max(0,muF)^3;
muFminus=max(0,-muF)^3;
lamFplus=max(0,lamF)^3;
lamFminus=max(0,-lamF)^3;
lamGplus=max(0,lamG)^3;
lamGminus=max(0,-lamG)^3;

% extract some other state-dependent variables
envec=instr.addvars;
bFpol=envec.bFpol;
bGpol=envec.bGpol;
kF=envec.kF;
levF=envec.levF;
kFpol = Kbar - kGpol;

% probabilities and states to compute expectation terms
prnext=obj.Exogenv.mtrans(exst,:);
a_next=obj.Exogenv.pts_all(:,1);
alower_next=obj.Exogenv.pts_all(:,2);

% projection evaluation
if nargin==7
    % Exactly 7 arguments were passed.
    % Means policies have already been computed but calcExpectations hasn't
	Pol_next = varargin{1};
elseif nargin==6
    % Nothing passed in ((conditional moments and EE errors in
    % simulation)
    Pol_next = obj.Vfct.evaluateAt(nextst)';
	if size(Pol_next,1)==1
        prnext=1;
	end
end

cF_next = Pol_next(:,1);
cG_next = Pol_next(:,2);
q_next = Pol_next(:,3);
kFpol_next = Pol_next(:,4);


% payoff of asset in worst state
if params.use_min_q
	min_q = min(q_next);
else
	min_q = q;
end


% FOCs for gatherer
MUg_next = betaG * (cG_next).^(-gamma);
exp_EEG_b = prnext*MUg_next;
exp_EEG_k = prnext*(MUg_next.*(alower_next*alpha.*kGpol.^(alpha-1) + q_next));

% FOCs for farmer
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

% Transitions
nst = length(a_next);
kFtrans = kFpol*ones(nst,1);
levFtrans = -R*bFpol./(q_next*kFpol);

V=cell(3,1);
J=[];
if mode==0 && nargout>1
	% Solution mode with analytic Jacobian requested
	
	% Prepare input
	cF_next_cell = num2cell(cF_next');
	cG_next_cell = num2cell(cG_next');
	q_next_cell = num2cell(q_next');
	kFpol_next_cell = num2cell(kFpol_next');
	prnext_cell = num2cell(prnext);
	
	a = instr.exstvec(1);
	alower = instr.exstvec(2);
	args = [ { R, q, cF, cG, kGpol, ...
				muFplus, muFminus, ...
				lamFplus, lamFminus, ...
				lamGplus, lamGminus, ...
				a, alower, kF, levF }, ...
				cF_next_cell, ...
				cG_next_cell, ...
				q_next_cell, ...
				kFpol_next_cell, ...
				prnext_cell];
			
	if params.use_min_q
		args = [args, {min_q}];
	end 
			
	% Compute jacobian
	fullJ = jacfun( args{:} );
    nsol = length(solvec);
    J = zeros(nsol);
	
	% Post-processing 1: apply chain rule to logged variables
	% Need: F'( y ) where y = log(x)
	% Have: F'( x ) = F'( exp(y) )
	% Chain rule: F'(y) = F'( exp(y) ) exp(y)
	startMult = 5;
	J(:,1:startMult) = fullJ(:,1:startMult) .* ...
		repmat( exp( solvec(1:startMult) )', size(fullJ,1), 1);
	
	% Need: F'(y) where y = max(0,muB)^3 or y = max(0,-muB)^3
	% Have: F'(x) = F'( max(0,muB)^3 ) or F'( max(0,-muB)^3 )
	% Chain rule: F'(y) = F'( max(0,muB)^3 ) * (3*muB^2) or 
	%					  F'( max(0,-muB)^3 ) * -(3*muB^2)
% 	if muF>0
% 		J = [J, fullJ(:,startMult) * (3 * muF^2) ];
% 	else
% 		J = [J, fullJ(:,startMult+1) * (-3 * muF^2)];
%     end
    
    for mult_idx = 1:nsol-startMult
        mult = solvec(startMult + mult_idx);
        if mult > 0
            J(:,startMult + mult_idx) = ...
                fullJ(:,startMult + 2*(mult_idx-1) + 1) .* repmat( 3*mult^2 ,nsol,1);
        else
            J(:,startMult + mult_idx) = ...
                fullJ(:,startMult + 2*(mult_idx-1) + 2) .* repmat( -3*mult^2 ,nsol,1);
        end
    end
    
elseif mode==1
    % Forecasting functions for time iteration
    Vnext=zeros(obj.Vfct.Nof,1);
    Vnext(1)=cF;
    Vnext(2)=cG;
    Vnext(3)=q;
    Vnext(4)=kFpol;
    V{1}=Vnext;
    % state transition
    V{2}=[kFtrans; levFtrans]';
    
elseif mode==2
    
    % Evaluation during simulation. Output conditional
    % variables
    
    % SDFs
	SDFG = MUg_next / MUg;
	SDFF = MUf_next / MUf;
    SDF.SDFG = SDFG;
    SDF.SDFF = SDFF;
    
    % Define returns
    retKF = (a_next + q_next)/q;
    retKG = (alower_next*alpha*kGpol^(alpha-1) + q_next)/q;
    expRF = prnext * retKF;
    expRG = prnext * retKG;
    
    condvars = struct('expRF',expRF, ...
        'expRG',expRG, ...
        'min_q',min_q);
    
    Wtrans.kF = kFtrans;
    Wtrans.levF = levFtrans;
    V = {condvars,Wtrans,SDF};
end

end
