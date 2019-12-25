clear;

close all;
respath='./';
outpath='./Results/';

% without and with adj. cost
if ~exist('econ', 'var'), econ='converge1'; end %econ='bench'; end
if ~exist('econ2', 'var'), econ2='theta875'; end
if ~exist('colors', 'var'), colors = {'b', 'r'}; end
%if ~exist('label', 'var'), labels = {'\alpha=0.7', '\alpha=0.9'}; end
if ~exist('label', 'var'), labels = {'\theta=0.9', '\theta=0.875'}; end
if ~exist('label', 'var'), titles = {'Capital F', 'Output'}; end

if exist('econ2', 'var')
    N_economy = 2;
else
    N_economy=1; % numbers economies to plot on same graph
end


if ~exist('resfile1', 'var'), resfile1=['res_',econ]; end
outpath=[outpath,'histcomp_',econ];
if N_economy == 2
    if ~exist('resfile2', 'var'), resfile2=['res_',econ2]; end

    outpath=[outpath,econ2];
end



load([respath,resfile1,'.mat']);
load([respath,'sim_',resfile1,'.mat']);
kFind=indexmap.get('kF');
Yind=indexmap.get('Y');

sim=simseries(:,[kFind,Yind]); 
sim(:,2)=sim(:,2)/mean(sim(:,2));
simseries_all{1} = sim; 
clear simseries;

if N_economy == 2
    load([respath,resfile2,'.mat']);
    load([respath,'sim_',resfile2,'.mat']);

    sim=simseries(:,[kFind,Yind]); 
    sim(:,2)=sim(:,2)/mean(sim(:,2));
    simseries_all{2} = sim;
    clear simseries;
end



%% histograms

datavec={{simseries_all{1}(:,1),simseries_all{2}(:,1)},{simseries_all{1}(:,2),simseries_all{2}(:,2)}};
edges{1}=linspace(quantile(simseries_all{2}(:,1),.001),quantile(simseries_all{2}(:,1),.99),18   );
edges{2}=linspace(quantile(simseries_all{2}(:,2),.001),quantile(simseries_all{2}(:,2),.99),18);
freqones=ones(size(datavec{1}{1}));

f=figure;
set(f,'PaperUnits','centimeters');
set(f,'PaperPosition',[2 10 20 10]);
hold on;
for i=1:2
    s=subplot(1,2,i);
    Helpers.makeHist2D(datavec{i},[],edges{i},titles{i},colors);
    xlims=[edges{i}(1)-abs(edges{i}(2)-edges{i}(1)),edges{i}(end)+abs(edges{i}(end)-edges{i}(end-1))];
    set(s,'Xlim',xlims);
    if i==1
        legend(labels,'Location','northwest');
    end
end

 
set(gcf,'PaperUnits','inches');
set(gcf,'PaperSize',[4*2 4*1]);
set(gcf,'PaperPosition',[0 0 4*2 4*1]);
set(gcf,'PaperPositionMode','manual');
%set(gcf,'PaperSize',[10*format(1), 3*format(2)]);
print('-dpdf','-r0',[outpath,'.pdf']);
print('-depsc',[outpath,'.eps']);





