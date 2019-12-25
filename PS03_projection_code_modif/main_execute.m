%========================================================================================
%   load model definition file
%========================================================================================

% open parallel pool
open_parpool(2); %had 10, but my pc only has 2

% path to file with experiment definition
exper_path='env_bench.mat';
% max iterations
maxit=200;
% max distance convergence
tol_max=1e-2;

% print mode
printmode=0;

% dampening
damp=0.5;

% revisiting of failed points
revisitFailed=true;

% load from disk
load(exper_path);
%load('res_2019_12_23_23_24.mat')


%========================================================================================
%   run policy iteration
%========================================================================================

gridSt=obj.Vfct.SSGrid.Pointmat; % use points from BaseGrid here
NDIM=obj.Vfct.SSGrid.Ndim;
exnpt=size(obj.Exogenv.pts_all,1);

% initialize
resmat=obj.Pfct.evaluateAt(gridSt)';
resmat_prev=resmat;

% split up matrix of points for better output
gr_points=cell(exnpt,1);
gr_index=cell(exnpt,2);
for i=1:exnpt
    grinlog=(gridSt(:,1)==i);
    grind=find(grinlog);
    gr_points{i}=gridSt(grinlog,:);
    gr_index{i,1}=grinlog;
    gr_index{i,2}=grind;
end

% value function
VF=obj.Vfct.evaluateAt(gridSt)';
VFnext=zeros(size(VF));
TF=obj.Tfct.evaluateAt(gridSt)';
TFnext=zeros(size(VF,1),obj.Tfct.Nof);

% control flags
iter=0;

disp(' ');
disp('Starting main loop ...');
disp(' ');
while 1
    % counter
    iter=iter+1;
    
    % ===========================================
    % loop over state space
    % ===========================================
    
    % matrix for failed points
    failedPoints=[];
    % transitions
    transmat=obj.Tfct.evaluateAt(gridSt)';
    
    failedPoints_trans_T=zeros(0,size(transmat,2));
    failedPoints_trans_I=[];
    failedPoints_trans_V=zeros(0,size(VF,2));
    
    % a rectangular grid to speed up VF interpolations solutions
    resmat_startiter = resmat; % for dampening
    % outer loop: all exogenous states
    for ei=1:exnpt
        tmp_grid=gr_points{ei};
        tmp_indlog=gr_index{ei,1};
        tmp_index=gr_index{ei,2};
        tmp_resmat=resmat(tmp_indlog,:);
        tmp_resmat_prev=resmat_prev(tmp_indlog,:);
        
        tmp_transmat=transmat(tmp_indlog,:);
        tmp_NPT=size(tmp_index,1); % nb SS pts for exog state ei
        
        % Evaluate value functions at transition points
        transpts=reshape([repmat(1:exnpt,tmp_NPT,1),tmp_transmat],tmp_NPT*exnpt,NDIM);
        Vtrans = obj.Vfct.evaluateAt(transpts)';
        
        % index matching for transitions
        refidx=kron((1:tmp_NPT)',ones(1,exnpt));
        refidx=refidx(:);
        
        
        disp(['State ',num2str(ei)]);
        
        [tmp_resmat_new,tmp_VF,tmp_TF,tmp_failed]=solvePointList(obj,tmp_grid,tmp_resmat,tmp_transmat,...
            refidx,Vtrans,tmp_resmat_prev,printmode,[]);
        failedPoints=[failedPoints; tmp_index(tmp_failed)];
        if revisitFailed
            failedPoints_trans_T=[failedPoints_trans_T; tmp_transmat(tmp_failed,:)];
            refidx_failed = ismember(refidx,find(tmp_failed));
            
            failedPoints_trans_I=[failedPoints_trans_I; tmp_index(refidx(refidx_failed))];
            failedPoints_trans_V=[failedPoints_trans_V; Vtrans(refidx_failed,:)];
        end
        
        resmat_prev(tmp_indlog,:)=tmp_resmat;
        resmat(tmp_indlog,:)=tmp_resmat_new;
        VFnext(tmp_indlog,:)=tmp_VF;
        TFnext(tmp_indlog,:)=tmp_TF;
    end
    
    if (revisitFailed && ~isempty(failedPoints))
        disp( '~~~~~~~~~~~~~~~~~~~');
        disp(['Revisiting failed points: ',num2str(length(failedPoints)),' add. points ...']);
        % try to solve at failed points
        [new_resmat,new_VF,new_TF,n_succ]=solvePointListFailed(obj,gridSt,failedPoints,resmat,...
            failedPoints_trans_T,failedPoints_trans_I,failedPoints_trans_V,...
            1,printmode,[]);
        resmat(failedPoints,:)=new_resmat;
        VFnext(failedPoints,:)=new_VF;
        TFnext(failedPoints,:)=new_TF;
    end
    
    
    % approximate functions for next iteration
    obj.Vfct=obj.Vfct.fitTo((1-damp)*VFnext+damp*VF);
    obj.Tfct=obj.Tfct.fitTo((1-damp)*TFnext+damp*TF);
    
    % convergence criterion (based on points in BaseGrid)
    val_range=1:3;
    VF_val = VF(:,val_range);
    VFnext_val=VFnext(:,val_range);
    [dist,wh]=max(abs(VF_val(:)-VFnext_val(:)));
    [mean_dist,col]=max(abs(mean(VF_val-VFnext_val)));
    [distT,whT]=max(abs(TF(:)-TFnext(:)));
    [wh_1,wh_2]=ind2sub(size(VFnext_val),wh);
    [whT_1,whT_2]=ind2sub(size(TFnext),whT);
    disp(['-- Iteration: ',num2str(iter),', max distance: ',num2str(dist),' in ',char(obj.V_names(val_range(1)-1+wh_2)), ...
        ' at point ',num2str(wh_1),': ',num2str(obj.Vfct.SSGrid.Pointmat(wh_1,:))]);
    disp(['-- Iteration: ',num2str(iter),', mean distance: ',num2str(mean_dist),' in ',char(obj.V_names(val_range(1)-1+col))]);
    disp(['-- Iteration: ',num2str(iter),', max T distance: ',num2str(distT),' in col ',num2str(whT_2), ...
        ' at point ',num2str(whT_1),': ',num2str(obj.Vfct.SSGrid.Pointmat(whT_1,:))]);
    disp(' ');
    if mean_dist<tol_max && dist<tol_max
        disp('Converged.');
        break;
    elseif iter>=maxit
        disp('Max.iter. exceeded.');
        break;
    end
    
    % update guess (based on points in BaseGrid)
    VF=VFnext;
    TF=TFnext;
end

% resulting policy functions
obj.Pfct=obj.Pfct.fitTo(resmat);



%========================================================================================
%   save result
%========================================================================================


% make date string
curr_date=clock;
date_str='';
for i=1:5
    date_str=[date_str,'_',num2str(curr_date(i))];
end

outname=['res',date_str,'.mat'];

save(outname,'obj','stv','failedPoints','dist','distT');




%========================================================================================
% infrastructure functions
%========================================================================================


% solution at point calling functions that contain model equations
function [fx,J,V]=solveAtPoint(obj,solvec,point,mode,trans,Polnext)
% transition
[nextst,outstr]=calcStateTransition(obj,point,solvec,0,trans);
% equations
[fx,J,V]=calcEquations(obj,point(1),nextst,solvec,outstr,mode,Polnext);
end

%-----------------------------------------------------------------------------------------

function [x,fx,exit,i]=tryOtherGuesses(fhand,gvec,options)
% list of guesses

%  Toggle Lagrange multipliers
gindex={6,6};
gvals={0.5,-0.5};

for i=1:length(gindex)
    newguess=gvec;
    newguess(gindex{i})=gvals{i};
    [x,fx,exit]=fsolve(fhand,newguess,options);
    if exit>0
        break;
    end
end
end

%-----------------------------------------------------------------------------------------

% solve model at list of points; called in main solution loop
function [resmat,VF,TF,failedPointsIndex]=solvePointList(obj,tmp_grid,tmp_resmat,tmp_transmat,...
    tmp_refidx,Vtrans,tmp_resmat_prev,print,logfname)

gr_npt=size(tmp_grid,1);
VF=zeros(gr_npt,obj.Vfct.Nof);
TF=zeros(gr_npt,obj.Tfct.Nof);
resmat=tmp_resmat;
failedPointsIndex=false(gr_npt,1);

% solver options
%             options=optimset('Display','off','TolX',1e-15,'TolFun',1e-12,...
%                 'MaxIter',60,'MaxFunEvals',60^2,'Algorithm','trust-region-dogleg');
options=optimset('Display','off','TolX',1e-15,'TolFun',1e-12,...
    'MaxIter',100,'MaxFunEvals',100^2,'FinDiffType','central',...
	'Jacobian','on');

f=Helpers.logopen(logfname);
fprintf(f,[repmat('.',1,gr_npt) '\n\n']);
Helpers.logclose(f);

%-----------------------------------------------------------------------------------------


% Use "for" for debugging. Should be "parfor" at runtime
parfor i=1:gr_npt
%for i=1:gr_npt
    x=[];fx=[];exfl=[];xret=[];exflret=[];count=[];
    ff=Helpers.logopen(logfname);
    
    point=tmp_grid(i,:);
    trans=tmp_transmat(i,:)';
    vtrans=Vtrans(tmp_refidx==i,:);
    
    succ='|';
    gvec=tmp_resmat(i,:)';
    
    % compute expectations
    transpts = [(1:obj.Exogenv.exnpt)',reshape(trans,obj.Exogenv.exnpt,obj.NSTEN)];
    Polnext = obj.Vfct.evaluateAt(transpts)';
	
    % prepare function handle
    fh_solveAtPoint=@(sol)solveAtPoint(obj,sol,point,0,trans,Polnext);
    
    % call to fsolve
    try
        [x,fx,exfl,~,J]=fsolve(fh_solveAtPoint,gvec,options);
		% for diagnostics only
		Jnum = jacobianest(fh_solveAtPoint,x);
		logmodulus = @(x)sign(x).*log(abs(x)+1);
		if max(max(abs(logmodulus(J)-logmodulus(Jnum))))>1e-4
			disp('Jacobian error');
		end
    catch E
        if print>0
            fprintf(ff,['Error at point ' ,num2str(i),':',num2str(point),'\n']);
        end
        rethrow(E);
    end
    
    if exfl<1 || sum(isnan(fx)) > 0
        % retry with different guess
        if print>1
            fprintf(ff,['Try other guesses at point ',num2str(i),':',num2str(point),'\n']);
        end
        try
            [xret,fx,exflret,count]=tryOtherGuesses(fh_solveAtPoint,gvec,options);
            
        catch E
            if print>0
                fprintf(ff,['Error (retry) at point ' ,num2str(i),':',num2str(point),'\n']);
            end
            rethrow(E);
        end
        if exflret>0
            xr=xret;
            succ=num2str(count);
        else
            if print>0
                fprintf(ff,['Return code: ',num2str(exfl),' at point ' ,num2str(i),':',num2str(point),'\n\n']);
            end
            xr=gvec;
            succ='x';
            failedPointsIndex(i)=true;
        end
    else
        xr=x;
    end
    
    % record result
    resmat(i,:)=xr;
    
    % also compute VF at solution
    [~,~,V]=solveAtPoint(obj,xr,point,1,trans,Polnext);
    
    VF(i,:)=V{1};
    TF(i,:)=V{2};
    fprintf(ff,'\b%s\n',succ);
    
    Helpers.logclose(ff);
end
end

%-----------------------------------------------------------------------------------------

% solve model at list of failed points
function [resmat_failed,VF_failed,TF_failed,n_succ]=solvePointListFailed(obj,grid,index_failed,resmat,...
    tmp_transmat, tmp_refidx,Vtrans, npasses,print,logfname)

gr_npt=size(grid,1);
resmat_failed=zeros(length(index_failed),size(resmat,2));
VF_failed=zeros(length(index_failed),obj.Vfct.Nof);
TF_failed=zeros(length(index_failed),obj.Tfct.Nof);

n_failed=length(index_failed);
idx_succ = zeros(length(index_failed));

% solver options
%             options=optimset('Display','off','TolX',1e-10,'TolFun',1e-12,...
%                 'MaxIter',100,'MaxFunEvals',20000,'Algorithm','trust-region-dogleg');
options=optimset('Display','off','TolX',1e-10,'TolFun',1e-12,...
    'MaxIter',100,'MaxFunEvals',100^2,'FinDiffType','central',...
	'Jacobian','on');


% keep track of failed indices
failed_idx_list=1:length(index_failed);

for p=1:npasses
    
    new_index_failed=false(size(index_failed));
    tmp_index_worked = setdiff(1:gr_npt,index_failed);
    tmp_grid_worked=grid(tmp_index_worked,:);
    
    f=Helpers.logopen(logfname);
    fprintf(f,[repmat('.',1,length(index_failed)) '\n\n']);
    Helpers.logclose(f);
    
    this_rmfailed=zeros(length(index_failed),size(resmat,2));
    this_VF_failed=zeros(length(index_failed),obj.Vfct.Nof);
    this_TF_failed=zeros(length(index_failed),obj.Tfct.Nof);
    % first loop over all failed points and find closest successful
    % point
	parfor i=1:length(index_failed)
    %for i=1:length(index_failed)
        ff=Helpers.logopen(logfname);
        exfl=0; x=0; fx=0;
        this_index=index_failed(i);
        point=grid(this_index,:);
        
        trans=tmp_transmat(i,:)';
        vtrans=Vtrans(tmp_refidx==this_index,:);
        
        % Find closest state
        dist=sum((tmp_grid_worked-repmat(point,length(tmp_grid_worked),1)).^2,2);
        [~,distidx]=sort(dist);
        ntr=1+2*p;
        gveclist = resmat(tmp_index_worked(distidx(1:ntr)),:);
        
        % compute expectations
        transpts = [(1:obj.Exogenv.exnpt)',reshape(trans,obj.Exogenv.exnpt,obj.NSTEN)];
		Polnext = obj.Vfct.evaluateAt(transpts)';;
    
		% prepare function handle
		fh_solveAtPoint=@(sol)solveAtPoint(obj,sol,point,0,trans,Polnext);
        for g=1:ntr
            gvec=gveclist(g,:)';
            % call to fsolve
            try
                [x,fx,exfl]=fsolve(fh_solveAtPoint,gvec,options);
                if exfl>0
                    break;
                end
            catch
                exfl=0;
            end
        end
        if exfl>0
            xr=x;
            succ='|';
            idx_succ(i)=1;
        else
            succ='x';
            new_index_failed(i)=true;
            xr=gveclist(1,:);
        end
        
        % record result
        this_rmfailed(i,:)=xr;
        
        % also compute VF at solution
        [~,~,V]=solveAtPoint(obj,xr,point,1,trans,Polnext);
        
        this_VF_failed(i,:)=V{1};
        this_TF_failed(i,:)=V{2};
        fprintf(ff,'\b%s\n',succ);
        
        Helpers.logclose(ff);
    end
    % write into return object
    resmat_failed(failed_idx_list,:)=this_rmfailed;
    VF_failed(failed_idx_list,:)=this_VF_failed;
    TF_failed(failed_idx_list,:)=this_TF_failed;
    failed_idx_list=failed_idx_list(new_index_failed);
    
    % update total resmat for guesses
    if sum(new_index_failed)==0 || sum(new_index_failed)==length(index_failed)
        break;
    end
    resmat(index_failed(~new_index_failed),:)=resmat_failed(~new_index_failed,:);
    index_failed=index_failed(new_index_failed);
    
end

n_succ = n_failed - length(index_failed);
end



