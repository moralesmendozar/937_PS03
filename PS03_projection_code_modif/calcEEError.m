function [errmat,solmat,condmat,Wshtrans,SDFmat]=calcEEError(obj,pointmat)
% function to compute Euler equation error at points in state
% space given by pointmat
nst=size(obj.Exogenv.pts_all,1);

errmat=zeros(size(pointmat,1),obj.Pfct.Nof);
solmat=zeros(size(errmat));
condmat=zeros(size(pointmat,1),obj.NCOND);
SDFmat=zeros(size(pointmat,1),2*nst);
Wshtrans=zeros(size(pointmat,1),2*nst);

evaluatePol = @(point)obj.Pfct.evaluateAt(point);
% Should be parfor. Use for when debugging only
for i=1:size(errmat,1)
    %            for i=1:size(errmat,1)
    point=pointmat(i,:);
    soltmp=evaluatePol(point);
    % transition
    [nextst,outstr]=calcStateTransition(obj,point,soltmp,0);
    % equations
    [fx,~,V]=calcEquations(obj,point(1),nextst,soltmp,outstr,2);
    R=exp(soltmp(1));
    q=exp(soltmp(2));
    bFpol=outstr.addvars.bFpol;
    Kbar=obj.Params.Kbar;
    normvec=[1/R,q,1/R,q,q,bFpol,Kbar,Kbar];
    condvars=V{1};
    kFtrans=V{2}.kF;
    levFtrans=V{2}.levF;
    SDFF=V{3}.SDFF;
    SDFG=V{3}.SDFG;
    errmat(i,:)=fx'./normvec;
    solmat(i,:)=soltmp';
    condmat(i,:)=Helpers.structToVec(condvars)';
    Wshtrans(i,:) = [kFtrans',levFtrans'];
    SDFmat(i,:) = [SDFF',SDFG'];
end

end
