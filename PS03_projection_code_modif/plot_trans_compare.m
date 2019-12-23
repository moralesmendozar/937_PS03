if ~exist('clear_flag', 'var'), clear_flag = 1; end

if usejava('desktop') && clear_flag
   clear; 
end

close all;
respath='./';
outpath='./Results/';

% without and with adj. cost
if ~exist('econ1', 'var'), econ1='bench'; end
if ~exist('econ2', 'var'), econ2='theta875'; end
if ~exist('colors', 'var'), colors = {'k-o', 'r-o'}; end
%if ~exist('label', 'var'), label = {'\alpha=0.7', '\alpha=0.9'}; end
if ~exist('label', 'var'), label = {'\theta=0.9', '\theta=0.875'}; end

if exist('econ3', 'var')
    N_economy = 3;
else
    N_economy=2; % numbers economies to plot on same graph
end


if ~exist('resfile','var')
    resfile=['res_comp_',econ1,econ2];
    if ~exist('resfile1', 'var'), resfile1=['res_',econ1]; end
    if ~exist('resfile2', 'var'), resfile2=['res_',econ2]; end
    outpath=[outpath,'GRcomp_',econ1,econ2];
    if N_economy == 3
        resfile = [resfile, econ3];
        if ~exist('resfile3', 'var'), resfile3 = ['res_',econ3]; end
        outpath=[outpath,econ3];
    end
    outpath=[outpath,'_'];
end
outfile=['GR_',resfile];
grayscale=0;
reportLevels=0;
batchMode=0;

load([respath,resfile1,'.mat']);
load([respath,'sim_',resfile1,'.mat']);
load([respath,'GR_',resfile1,'.mat']);

simseries_mean_econ{1} = simseries_mean{3}; % financial recession 
clear simseries_mean;

load([respath,resfile2,'.mat']);
load([respath,'sim_',resfile2,'.mat']);
load([respath,'GR_',resfile2,'.mat']);

simseries_mean_econ{2} = simseries_mean{3};
clear simseries_mean;

if N_economy == 3
    load([respath,resfile3,'.mat']);
    load([respath,'sim_',resfile3,'.mat']);
    load([respath,'GR_',resfile3,'.mat']);

    simseries_mean_econ{3} = simseries_mean{3};
    clear simseries_mean;
end

close all;

%% IRFS

% file names for graphs (set to empty for no printing)
printfiles={[outpath,'IRF1'],[outpath,'IRF2']};        
% which variables
brsel1=[indexmap.get('a'),indexmap.get('Y'),indexmap.get('kF'),...    
         indexmap.get('cF'),indexmap.get('cG'),indexmap.get('wFsh')];
levelind1=zeros(size(brsel1));     
brsel2=[indexmap.get('R'),indexmap.get('q'),...
        indexmap.get('bFpol'),indexmap.get('bind_muF')];
levelind2=zeros(size(brsel2));     

titles1={'TFP','Output','Cap F','Cons F','Cons G','Wealth F'}; 
titles2={'R','q','Debt F','\muF'};

tvec=0:NT_sim-1;    
    
brsel_all=[brsel1,brsel2];
levelind_all=[levelind1, levelind2];
nvar=length(brsel_all);   
brseries_gr=zeros(N_economy, NT_sim, nvar);
for s=1:N_economy
    for v=1:nvar
        thisv=brsel_all(v);
        if levelind_all(v) 
            brseries_gr(s,:,v) = simseries_mean_econ{s}(:,thisv);
        else
            brseries_gr(s,:,v) = 100*(simseries_mean_econ{s}(:,thisv) ./ repmat(simseries_mean_econ{s}(1,thisv),NT_sim,1)-1);
        end
    end
end


if usejava('desktop')
    Helpers.makeImpRes(brseries_gr(:,:,1:6),tvec,titles1,colors,[2,3],label,printfiles{1});   
    Helpers.makeImpRes(brseries_gr(:,:,7:10),tvec,titles2,colors,[2,2],label,printfiles{2});    
end



