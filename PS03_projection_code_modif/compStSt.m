function stvals=compStSt(params,print)
% unpack relevant parameters
betaF=params.betaF;
betaG=params.betaG;
alpha=params.alpha;
a=params.a;
alower=params.alower;
Kbar=params.Kbar;
theta=params.theta;

%muF = betaG - betaF;
muF = 1 - betaF/betaG;
% solution for KG
%kGunc = (alower/a * betaG/betaF * alpha*(1- (1-theta)*betaF - theta*betaG)/(1-betaG) )^(1/(1-alpha));
kGunc = ( (alpha * alower * betaG * (1 - betaF - theta*muF ) )/...
    ((1-betaG)*a*(betaF + theta*muF)) )^(1/(1-alpha));
lamG = 0;
lamF = 0;
if kGunc <= 0
    kG = 0;
    lamG = -kGunc;
elseif kGunc >= Kbar
    kG = Kbar;
    lamF = kGunc - Kbar;
else
    kG = kGunc;
end
kF = Kbar - kG;
R = 1/betaG;
%q = ( (betaF)*a )/(1-betaF-theta*muF);
q = ( (betaF+theta*muF)*a )/(1-betaF-theta*muF);
%bF = -theta*q*kF;
bF = -theta*(q+a)*kF;
bG = -bF;
yF = a*kF;
yG = alower*kG^alpha;
Y = yF + yG;
wF = (a+q)*kF + R*bF;
wFsh = (q*kF + R*bF)/(q*Kbar);
cF = wF - q*kF - bF;
wG = yG + q*kG + R*bG;
cG = wG - q*kG - bG;
levF= -R*bF/(q*kF);


if print
    % print steady state values
    disp(' ');
    disp('Deterministic steady state');
    disp('--- Prices ---');
    disp(['R: ',num2str(R)]);
    disp(['q: ',num2str(q)]);
    disp(['muF: ',num2str(muF)]);
    disp(['lamF: ',num2str(lamF)]);
    disp(['lamG: ',num2str(lamG)]);
    disp('--- Output ---');
    disp(['Y: ',num2str(Y)]);
    disp(['yF: ',num2str(yF)]);
    disp(['yG: ',num2str(yG)]);
    disp('--- Wealth ---');
    disp(['kF: ',num2str(kF)]);
    disp(['wF: ',num2str(wF)]);
    disp(['wFsh: ',num2str(wFsh)]);
    disp(['levF: ',num2str(levF)]);
    disp(['wG: ',num2str(wG)]);
    disp(['bF: ',num2str(bF)]);
    disp(['bG: ',num2str(bG)]);
    disp('--- Consumption ---');
    disp(['cF: ',num2str(cF)]);
    disp(['cG: ',num2str(cG)]);
end

Sol=struct('R',R,...
    'q',q,...
    'cF',cF,...
    'cG',cG,...
    'kGpol',kG,...
    'muF',sign(muF)*abs(muF)^(1/3),...
    'lamF',lamF^(1/3),...
    'lamG',lamG^(1/3));

Soldomain=struct('R',[0,Inf],...
	'q',[0,Inf],...
	'cF',[0,Inf],...
	'cG',[0,Inf],...
	'kGpol',[0,Inf],...
	'muF',[-Inf,Inf],...
    'lamF',[-Inf,Inf],...
    'lamG',[-Inf,Inf]);


V=struct('cF',cF,...
    'cG',cG,...
    'q',q,...
    'kFpol',kF);

Add=struct('bFpol',bF,...
    'bGpol',bG,...
    'kGpol',kG,...
    'kF',kF,...
    'levF',levF,...
    'Y',Y,...
    'wG',wG,...
    'wF',wF);


State=struct('kF',kF,...
    'levF',levF);

stvals=struct('Sol',Sol,...
	'Soldomain',Soldomain,...
    'V',V,...
    'Add',Add,...
    'State',State);

end
