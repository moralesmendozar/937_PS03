clear; 
close all;

%--------------------------------------------------------------------------
% simulation setup
%--------------------------------------------------------------------------

% file with model
respath='./';
outpath='./Results/';
%resfile_list={'res_theta875'}; 
resfile_list={'res_01'}; 

for f=1:length(resfile_list)
    
    resfile=resfile_list{f};
    
    load([respath,resfile,'.mat']);
    
    % Initial Economy Config
    varlist={'simseries','statevec','indexmap','varnames'};
    load(['sim_',resfile],varlist{:});
    
    % set starting point
    start_ini=3;
    start_shock=[0,2,1]; 
    statevec=statevec(2:end);
    startvals=mean(simseries(statevec==start_ini,:));
    N_vars=length(startvals);
        
    % number of periods and burn-in
    N_shock=length(start_shock);
    N_runs=5000;
    NT_sim=25;
    NT_ini=0;
    
    % cluster endogenous states
    envarind=3:4;
    maxclust=10;
    enstatemat=simseries(statevec==start_ini,envarind);
    cindex=clusterdata(enstatemat,'criterion','distance','maxclust',maxclust,'linkage','weighted');
    sttot=size(enstatemat,1);
    startptmat=[];
    for c=1:maxclust
        cfrac=sum(cindex==c)/sttot;
        thismean=mean(enstatemat(cindex==c,:),1);
        disp([num2str(c),': ',num2str(thismean)]);
        thisc=repmat(thismean,floor(N_runs*cfrac),1);
        startptmat=[startptmat; thisc];
    end
    if size(startptmat,1)<N_runs
        thismean=mean(enstatemat);
        startptmat=[startptmat; repmat(thismean,N_runs-size(startptmat,1),1)];
    end
    
    
    % report levels or grwoth rates for output variables
    reportLevels=1;
    
    % compute Euler equation error?
    compEEErr=1;
    % make graphs grayscale
    grayscale=0;
    
    % Compute term premium (slow!)
    term_premium=0;
    
    % output table file
    outfile=['GR_',resfile];
    
    varnames_store = varnames;
    
    simseries_median = cell(N_shock,1);
    simseries_mean = cell(N_shock,1);
    simseries_std = cell(N_shock,1);
    
    simseries_diff_median = cell(N_shock,1);
    simseries_diff_mean = cell(N_shock,1);
    simseries_diff_std = cell(N_shock,1);
    
    open_parpool(2);
        
    for s=1:N_shock
        
        disp(['Shock ',num2str(s),' of ',num2str(N_shock)]);
        tens_simseries = zeros(NT_sim+1,N_vars,N_runs);
        
        % compute entry of random number matrix that sets first state
        % deterministically to start_shock
        if start_shock(s)>0
            transprob=cumsum(obj.Exogenv.mtrans(start_ini,:));
            shock_prob=transprob(start_shock(s));
            if start_shock(s)>1
                shock_prob_minus=transprob(start_shock(s)-1);
            else
                shock_prob_minus=0;
            end
            rvar_next=(shock_prob+shock_prob_minus)/2;
        end
        
        % Create shock matrix
        rng(1);
        shmatfull = lhsdesign(N_runs,NT_sim);
        
        SDFmat=zeros(NT_sim,obj.Exogenv.exnpt);
        
        fprintf([repmat('.',1,100) '\n\n']);
        
        parfor n=1:N_runs
            %--------------------------------------------------------------------------
            % start simulation
            %--------------------------------------------------------------------------
            %fprintf('Run %d - Start \n',n);
            % simulate
            shmat = shmatfull(n,:)';
            if start_shock(s)>0
                shmat(1)=rvar_next;
            end
            
            startpt=struct;
            startpt.kF=startptmat(n,1);
            startpt.levF=startptmat(n,2);
            startpt=orderfields(startpt,obj.En_names);
            startpt_vec=Helpers.structToVec(startpt)';
            startpt_vec=[start_ini,startpt_vec];
            
            [simseries,varnames,~,~,~]=simulate(obj,NT_sim,NT_ini,startpt_vec,compEEErr,shmat);
            simseries_orig=simseries;
            varnames_orig=varnames;
            statevec = simseries(:,1);
            %fprintf('Run %d - After simulation \n',n);
            
            [simseries, varnames] = computeSimulationMoments(simseries,varnames);
            
            nvars = length(varnames);
            %fprintf('Run %d - After computation \n',n);
            %         disp(size(startvals))
            %         disp(size(simseries))
            tens_simseries(:,:,n) = [startvals; simseries];
            if mod(n,N_runs/100)==0
                fprintf('\b|\n');
            end
            
        end
        fprintf('\n');
        varnames = varnames_store;
        nvars = length(varnames);
        
        % make HashMap with mapping of names to indices
        indexmap=java.util.HashMap;
        for i=1:nvars
            indexmap.put(varnames{i},i);
        end

        
        simseries_median{s} = median(tens_simseries,3);
        simseries_mean{s} = mean(tens_simseries,3);
        simseries_std{s} = std(tens_simseries,[],3);
        
        if start_shock(s) > 0
            % If actual shock, difference and save
            tens_simseries_diff = tens_simseries - tens_simseries_0;
            simseries_diff_median{s} = median(tens_simseries_diff,3);
            simseries_diff_mean{s} = mean(tens_simseries_diff,3);
            simseries_diff_std{s} = std(tens_simseries_diff,[],3);
        else
            % If no shock, store for later differencing
            tens_simseries_0 = tens_simseries;
        end

        
    end
    
    NT_sim=NT_sim+1;
    save(outfile,'simseries_mean','simseries_median','simseries_std', ...
        'simseries_diff_mean', 'simseries_diff_median', 'simseries_diff_std', ...
        'indexmap','NT_sim','N_shock');
    
end
