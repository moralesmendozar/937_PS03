function [simseries,varnames,errmat,Wshtrans,SDFmat]=simulate(obj,NT,NTini,inistvec,simerror,shmat)
if length(inistvec)~=obj.Vfct.SSGrid.Ndim
    error('inistvec must be vector of length SSGrid.Ndim');
end

NTtot=NT+NTini;
simseries=zeros(NTtot,1+obj.NSTEX+obj.NSTEN+obj.NSOL+obj.NV+obj.NADD);

% if shock matrix wasn't passed in, create it
preset_path = false;
if nargin<6
    rng(10,'twister');
    shmat=lhsdesign(NTtot,1,'criterion','correlation'); % realizations of shocks for Markov state
else
    preset_path=isinteger(shmat);
end

point=inistvec;

pointmat=zeros(NTtot,length(point));

for t=1:NTtot
    pointmat(t,:)=point;
    exst=point(1);
    
    % next period's exog. state
    if preset_path
        exnext=shmat(t);
    else
        transprob=cumsum(obj.Exogenv.mtrans(exst,:));
        exnext=find(transprob-shmat(t)>0,1,'first');
    end
    
    % transition to next period
    solvec=obj.Pfct.evaluateAt(point);
    valvec=obj.Vfct.evaluateAt(point)';
    [nextst,outstr]=calcStateTransition(obj,point,solvec,exnext);
    
    addvec=Helpers.structToVec(outstr.addvars)';
    % write different categories of variables in one row
    simnext=[point(1),outstr.exstvec',point(2:end),solvec',valvec,addvec];
    if length(simnext)~=size(simseries,2)
        disp('problem');
    end
    simseries(t,:)=simnext;
    point=nextst;
end

simseries=simseries(NTini+1:end,:);
varnames=[{'exst'}, obj.Ex_names, obj.En_names, obj.Sol_names, obj.V_names, obj.Add_names];

errmat=[];
Wshtrans=[];%zeros(size(obj.Exogenv.mtrans(simseries(:,1),:)));
SDFmat=[];%zeros(size(obj.Exogenv.mtrans(simseries(:,1),:)));
if simerror 
    [errmat,~,condmat,Wshtrans,SDFmat]=calcEEError(obj,pointmat);
    errmat=errmat(NTini+1:end,:);
    condmat=condmat(NTini+1:end,:);
    Wshtrans=Wshtrans(NTini+1:end,:);
    SDFmat=SDFmat(NTini+1:end,:);
    simseries=[simseries,condmat];
    varnames=[varnames,obj.Cond_names];
else
    simseries=[simseries,zeros(NT,obj.NCOND)];
    varnames=[varnames,strcat(obj.Cond_names,'_nan')];
end
end
