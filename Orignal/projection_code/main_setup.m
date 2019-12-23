function [stv,obj] = main_setup(experfile,expername,guess_path,useJac) 

if nargin<4
	useJac=true;
	if nargin<3
		% if using an existing obect file as guess,
		% set path here, otherwise set to empty
		%guess_path='res_phi03.mat';
		guess_path=[];
		if nargin<2
			% Define default experiment
			expername = 'bench';
			if nargin<1
				% Define default experdef file
				experfile = 'experdef';
			end
		end
	end
end

tmp = str2func(experfile);
allexpers = tmp();
params = allexpers.(expername).params;
grids = allexpers.(expername).grids;

%========================================================================================
%   compute steady state
%========================================================================================

stv=compStSt(params,1);

solguess = stv.Sol;
for varname = fieldnames(solguess)'
	varname = varname{:};
	if sum( stv.Soldomain.(varname) == [0,Inf] ) == 2
		solguess.(varname) = log(stv.Sol.(varname));
	end
end

% Modify initial guess to differ from deterministic steady state
solguess.muF = 0.5;
solguess.lamF = -0.5;
solguess.lamG = -0.5;

solguessvec=Helpers.structToVec(solguess);
Vguessvec=Helpers.structToVec(stv.V);
V_names=fieldnames(stv.V);


% ===================================================
% Parameters of stochastic model
% ===================================================

% for growth rate of income

params.mu_a=-.5*params.sig_a^2/(1-params.rho_a^2);
Skew_a=0;
[aprob,avec] = Helpers.rouwen(params.rho_a,params.mu_a,params.sig_a, ...
								Skew_a,params.N_a);
aprob=aprob';

mpts_perm=[exp(avec),exp(avec)*params.alower];
exnpt=size(mpts_perm,1);
mpts_all=[mpts_perm, (1:exnpt)'];
     

%========================================================================================
%   create model structure with initial guess
%========================================================================================

% exogenous state variables
exogenv=struct;
exogenv.exnames={'a','alower'};
exogenv.exnpt=exnpt;
exogenv.pts_all=mpts_all;
exogenv.mtrans=aprob;

% endogenous model objects
ennames=fieldnames(grids);                              
solnames=fieldnames(stv.Sol);
solbase=solguessvec;
addnames=fieldnames(stv.Add);
condnames = {'expRF', 'expRG', 'min_q'};                        

% make grid object
unigrids= [ {1:exogenv.exnpt}, struct2cell(grids)' ];
basegrid=grid.TensorGrid(unigrids);

if isempty(guess_path)
    solmatguess=repmat(solguessvec',[basegrid.Npt,1]);
    Vmatguess=repmat(Vguessvec',[basegrid.Npt,1]);
    transguess=kron(basegrid.Pointmat(:,2:end), ones(1,exogenv.exnpt));

	%Vmatguess(:,4) = basegrid.Pointmat(:,2);
else
    guessobj=load(guess_path,'obj');
    solmatguess=guessobj.obj.Pfct.evaluateAt(basegrid.Pointmat)';
    Vmatguess=guessobj.obj.Vfct.evaluateAt(basegrid.Pointmat)';
    transguess=guessobj.obj.Tfct.evaluateAt(basegrid.Pointmat)';
end
            
% build approximating functions
% create guess for solution variable functions
Pf=grid.LinearInterpFunction(basegrid,solmatguess);
% create guess for next-period functions
Vf=grid.LinearInterpFunction(basegrid,Vmatguess);
% and for state transition function
Tf=grid.LinearInterpFunction(basegrid,transguess);
           
% create new model environment structure
obj=struct;
obj.Pfct=Pf;
obj.Vfct=Vf;
obj.NV=Vf.Nof;
obj.V_names=reshape(V_names,1,obj.NV);
obj.Tfct=Tf;
obj.NSOL=length(solnames);
obj.Sol_names=reshape(solnames,1,obj.NSOL);
if obj.NSOL~=Pf.Nof
    error('Inconsistent number of solution variables.');
end
obj.Sol_baseguess=solbase;
obj.NSTEN=length(ennames);
obj.En_names=reshape(ennames,1,obj.NSTEN);
if obj.NSTEN~=Vf.SSGrid.Ndim-1
    error('Inconsistent number of endogenous state variables.');
end
obj.NADD=length(addnames);
obj.Add_names=reshape(addnames,1,obj.NADD);
obj.Params=params;
obj.Exogenv=exogenv;
obj.NSTEX=length(exogenv.exnames);
obj.Ex_names=reshape(exogenv.exnames,1,obj.NSTEX);
obj.NCOND=length(condnames);
obj.Cond_names=condnames;

% Analytic jacobian
if useJac
	obj.Jac = constructJacobian(obj);
else
	obj.Jac = [];
end

disp('-------------------------------------');
disp('Bounds of state space:');
disp(num2str(obj.Vfct.SSGrid.StateBounds));

% ===========================================
% save initial model definition
% ===========================================

outfname=['env_',expername,'.mat'];
save(outfname,'stv','obj');              




