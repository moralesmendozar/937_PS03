clear; 
close all;

respath='./';
outpath='./Results/';
if ~exist('resfile','var')
    %resfile='res_phi04';
    resfile='res_01';
end
outfile=['GR_',resfile];
grayscale=0;
reportLevels=0;
batchMode=0;
relative_irf=1;

load([respath,resfile,'.mat']);
load([respath,'sim_',resfile,'.mat']);
load([respath,'GR_',resfile,'.mat']);

close all;

% Make  Graphs
outpath=[outpath,outfile,'_'];
if relative_irf, outpath = [outpath, 'diff_']; end
tvec=0:NT_sim-1;

% file names for graphs (set to empty for no printing)
printfiles={[outpath,'IRF1'],[outpath,'IRF2']};        
% which variables
brsel1=[indexmap.get('a'),indexmap.get('Y'),indexmap.get('kF'),...
         indexmap.get('cF'),indexmap.get('cG'),indexmap.get('wFsh')];
brsel2=[indexmap.get('R'),indexmap.get('q'),...
        indexmap.get('bFpol'),indexmap.get('bind_muF')];

% How many shocks to plot (one less for relative)
if relative_irf
    N_shock_plot = N_shock - 1;
else
    N_shock_plot = N_shock;
end

brsel_all=[brsel1,brsel2];
nvar=length(brsel_all);   
brseries_gr=zeros(N_shock_plot, NT_sim, nvar);

for s=1:N_shock_plot
    if relative_irf
        if reportLevels==1
            brseries_gr(s,:,:) = simseries_diff_mean{s+1}(:,brsel_all);
        else
            brseries_gr(s,:,:) = 100*(simseries_diff_mean{s+1}(:,brsel_all) ./ repmat(simseries_mean{s+1}(1,brsel_all),NT_sim,1));
        end
    else
        if reportLevels==1
            brseries_gr(s,:,:) = simseries_mean{s}(:,brsel_all);
        else
            brseries_gr(s,:,:) = 100*(simseries_mean{s}(:,brsel_all) ./ repmat(simseries_mean{s}(1,brsel_all),NT_sim,1)-1);
        end
    end
end

colors={'b-o','r-o'};
if N_shock_plot==3
   colors = ['k-o',colors]; 
end
titles1={'TFP','Output','Cap F','Cons F','Cons G','Wealth F'}; 
titles2={'R','q','Debt F','\muF'};
     
if usejava('desktop')
    Helpers.makeImpRes(brseries_gr(:,:,1:6),tvec,titles1,colors,[2,3],[],printfiles{1});   
    Helpers.makeImpRes(brseries_gr(:,:,7:10),tvec,titles2,colors,[2,2],[],printfiles{2});
end


if batchMode
    close all;
    clear;
end
