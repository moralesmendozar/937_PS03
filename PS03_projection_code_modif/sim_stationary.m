clear; 
close all;

%--------------------------------------------------------------------------
% simulation setup
%--------------------------------------------------------------------------

% file with model
respath='./';
%resfile='res_alpha90'; 
resfile='res_converge1'; 
load([respath,resfile,'.mat']);

% Update params
augmentParams=0;
expdef='experdef_20181029.m';


% number of periods and burn-in
NT_sim=10000;
NT_ini=500;

% compute Euler equation error?
compEEErr=1;

% Force the creation of a sim_res file
force_output = 1;

% output table file
output_dir = './Results/';
if ~exist(output_dir, 'dir')
    mkdir(output_dir)
end

outstats_exog=[output_dir,'statsexog_',resfile,'.xls'];
errstats=[output_dir,'errstats_',resfile,'.xls'];
expdata=0;
outdata=[output_dir,'series_',resfile,'.csv'];
       
%--------------------------------------------------------------------------
% start simulation
%--------------------------------------------------------------------------

% set starting point
start_ex=3;
startpt=struct;
startpt.kF=stv.State.kF;
startpt.levF=stv.State.levF;
startpt=orderfields(startpt,obj.En_names);
startpt_vec=Helpers.structToVec(startpt)';
startpt_vec=[start_ex,startpt_vec];

% simulate
[simseries,varnames,errmat,Wshtrans,SDFmat]=simulate(obj,NT_sim,NT_ini,startpt_vec,compEEErr);
simseries_orig=simseries;
varnames_orig=varnames;
statevec = simseries(:,1);


[simseries, varnames] = computeSimulationMoments(simseries,varnames);
nvars = length(varnames);

% Create table object for easier access
simtable=array2table(simseries);
[~,ia,~]=unique(varnames);
simtable=simtable(:,ia);
simtable.Properties.VariableNames=varnames(ia);
dispnames=varnames(ia);

% make HashMap with mapping of names to indices
indexmap=java.util.HashMap;
for i=1:nvars
    indexmap.put(varnames{i},i);
end

% Check transition function errors
if compEEErr
    idx = sub2ind([NT_sim,obj.Exogenv.exnpt],(1:NT_sim)',[statevec(2:end);1]);
    idx=idx(1:end-1);
    kFtrans=Wshtrans(:,1:obj.Exogenv.exnpt);
    kF_err=simseries(2:end,indexmap.get('kF')) - kFtrans(idx); % column index: 6
    errmat = [errmat, [kF_err;0]];
    levFtrans=Wshtrans(:,obj.Exogenv.exnpt+1:end);
    levF_err=simseries(2:end,indexmap.get('levF')) - levFtrans(idx);
    errmat = [errmat, [levF_err;0]];
end

%--------------------------------------------------------------------------
% calculate stats
%--------------------------------------------------------------------------
varst=zeros(length(startpt_vec)-1,1);
for i=1:length(startpt_vec)-1
    varst(i)=indexmap.get(obj.En_names{i});
end
            
% state variable means in stationary distribution
stvstat=mean(simseries(:,varst));

% calculate business cycle stats
statsout_exog=cell(4,1);
statsout_endog=cell(4,1);

% first for all periods, then separately for good and bad states
smpsel_exog={true(NT_sim-1,1), simseries(:,1) >= obj.Params.a, simseries(:,1) < obj.Params.a};


gdp_idx = indexmap.get('Y');                     
for j=1:numel(smpsel_exog)
    % subsample
    simtmp=simseries(smpsel_exog{j},:);
    statstmp=zeros(nvars,7);
    statstmp(:,1)=nanmean(simtmp)';
    statstmp(:,2)=nanstd(simtmp)';
    % contemp and first-order autocorrelations
    autocorrm=corrcoef([simtmp(2:end,:),simtmp(1:end-1,:)]);
    conm=autocorrm(1:nvars,1:nvars);
    lagm=autocorrm(nvars+1:end,1:nvars);
    % corr with shocks
    statstmp(:,3:4)=[conm(:,1),lagm(1,:)'];
    % corr with Y
    statstmp(:,5:6)=[conm(:,gdp_idx),lagm(gdp_idx,:)'];
    % vector with fo autocorr
    statstmp(:,7)=diag(lagm);
    statsout_exog{j}=statstmp;
    
end


%--------------------------------------------------------------------------
% output
%--------------------------------------------------------------------------

% overview output for eyeball check against analytic st.st. values
% make one big structure with steady-state values
stvbig=Helpers.combineStructs({stv.Sol,stv.State,stv.Add});

% output table
% make index vector
[displist,dispnames]=Helpers.makeListFromNames(indexmap,dispnames);
ndvars=length(displist);

disp(' ');
disp('Simulation steady state');

% overview output 
fprintf('Frequency (exog subsamples): ');
for j=1:numel(smpsel_exog)
    % select vars
    tabout_exog{j}=statsout_exog{j}(displist,:);
    fprintf('%f\t',sum(smpsel_exog{j}));
end
fprintf('\n');
disp('-------------');

for s=1:ndvars
    if isfield(stvbig,dispnames{s})
        ststval=stvbig.(dispnames{s});
    else
        ststval=0;
    end
    if numel(dispnames{s}) > 7
        fprintf('%d\t%4s\t\t\t\t%f |',displist(s),dispnames{s},ststval);
    else
        fprintf('%d\t%4s\t\t\t\t\t%f |',displist(s),dispnames{s},ststval);
    end
%     disp('Exog subsamples')
    for j=1:numel(smpsel_exog)
        fprintf('\t%f, %f |',tabout_exog{j}(s,1),tabout_exog{j}(s,2));
    end
    fprintf('\n');    
end

if compEEErr
    avg_err=mean(abs(errmat))';
    med_err=median(abs(errmat))';
    p75_err=prctile(abs(errmat),75)';
    p95_err=prctile(abs(errmat),95)';
    p99_err=prctile(abs(errmat),99)';
    p995_err=prctile(abs(errmat),99.5)';
    max_err=max(abs(errmat))';
    errtab=table(avg_err,med_err,p75_err,p95_err,p99_err,p995_err,max_err);
    errarr=table2array(errtab);
    disp(' ');
    disp('-----------------------------------------------');
    disp('Average and maximum Euler equation error');
    fprintf('Equ.no.\t\tAvg.\t\tMed.\t\tp75\t\t\tp95\t\t\tp99\t\t\tp99.5\t\tMax.\n');
    for s=1:length(avg_err)
        fprintf('%d\t\t\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n',s,errarr(s,1),errarr(s,2),errarr(s,3), ...
            errarr(s,4),errarr(s,5),errarr(s,6),errarr(s,7));
    end
        
end

% check grid bounds
state_range=3:4;
min_vec=min(simseries(:,state_range));
max_vec=max(simseries(:,state_range));
disp('State bounds:');
disp(obj.Pfct.SSGrid.StateBounds(:,2:end));
disp('Simulation mins:');
disp(min_vec);
disp('Simulation max:');
disp(max_vec);


% write to file
values=struct2cell(obj.Params);
paramout=cell2table(values,'RowNames',fieldnames(obj.Params));
colnames={'mean','std','corrG','corrG_1','corrY','corrY_1','AC'};
for j=1:numel(smpsel_exog)
    tableout_exog=array2table(tabout_exog{j},'RowNames',dispnames,'VariableNames',colnames);
    writetable(tableout_exog,outstats_exog,'WriteRowNames',1,'FileType','spreadsheet','Sheet',j);
end
writetable(paramout,outstats_exog,'WriteRowNames',1,'FileType','spreadsheet','Sheet','params');
writetable(errtab,errstats,'FileType','spreadsheet');
if force_output
    params=obj.Params;
    disp(['Saving simulation data to .mat file: ',['sim_',resfile,'.mat']]);
    save(['sim_',resfile,'.mat'],'simseries','displist','dispnames','errmat','tabout_exog','outstats_exog', ...
                'errstats','errtab','indexmap','NT_ini','NT_sim','smpsel_exog','statevec','statsout_exog','varnames','params');
end    

if expdata
    disp(' ');
    disp('Exporting simseries...');
    Helpers.tableExport(outdata,varnames,simseries);
end

% save model file with stationary state values
save([respath,resfile,'.mat'],'stvstat','-append');

