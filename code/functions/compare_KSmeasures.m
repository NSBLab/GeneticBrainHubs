% compare ST and GT energy distributions
% 1 - is ST
% 3 - GT
DATA = cell(5,1); 
DATA{1} = load('Group_HCPparc20dens_10k_space1_deg-avg1_gene0_voronoi_mult2_energy_sp.mat'); 
DATA{2} = load('Group_HCPparc20dens_10k_space1_sptl1_gene1_voronoi_mult2_energy_sp.mat'); 
DATA{3} = load('Group_HCPparc20dens_10k_space0_deg-avg1_gene1_voronoi_mult2_energy_sp.mat'); 
DATA{4} = load('Group_HCPparc20dens_10k_space1_sptl1_gene0_voronoi_mult2_energy_sp.mat'); 
DATA{5} = load('Group_HCPparc20dens_10k_space0_sptl1_gene1_voronoi_mult2_energy_sp.mat'); 

% label each model in the same order
% S-stands for space; T - topology; G - gene; 
mtype = {'ST', 'SG', 'TG', 'S', 'G'}; 

ENERGY = cell(1,length(mtype)); 

for m=1:length(mtype)

    ENERGY{m} = DATA{m}.E;

end

[S, INDselected] = plotMODviolin(ENERGY, 100, mtype, 'lowest');

P.('ST') = DATA{1}.K(INDselected{1},:); 
P.('TG') = DATA{3}.K(INDselected{3},:); 

figure; 
subplot(1,4,1); histogram(P.ST(:,1)); hold on; histogram(P.TG(:,1)); title('Degree'); xlabel('KS'); xlim([0 0.25])
subplot(1,4,2); histogram(P.ST(:,2)); hold on; histogram(P.TG(:,2)); title('Clust'); xlabel('KS'); xlim([0 0.25])
subplot(1,4,3); histogram(P.ST(:,3)); hold on; histogram(P.TG(:,3)); title('Betweenness'); xlabel('KS'); xlim([0 0.25])
subplot(1,4,4); histogram(P.ST(:,4)); hold on; histogram(P.TG(:,4)); title('Length'); xlabel('KS'); xlim([0 0.25])
legend({'ST', 'TG'})