function Figure5()

hemi = 'left'; %if using degCorr, select only 'left'
parcellation = 'HCP';
optimiseWhat = 'energy';

% load all data data 
% sptl doesn't have any topological iinformation
DATA = cell(5,1); 
DATA{1} = load('Group_HCPparc20dens_10k_space1_deg-avg1_gene0_voronoi_mult2_energy_sp.mat'); 
DATA{2} = load('Group_HCPparc20dens_10k_space1_sptl1_gene1_voronoi_mult2_energy_sp.mat'); 
DATA{3} = load('Group_HCPparc20dens_10k_space0_deg-avg1_gene1_voronoi_mult2_energy_sp.mat'); 
DATA{4} = load('Group_HCPparc20dens_10k_space1_sptl1_gene0_voronoi_mult2_energy_sp.mat'); 
DATA{5} = load('Group_HCPparc20dens_10k_space0_sptl1_gene1_voronoi_mult2_energy_sp.mat'); 

% label each model in the same order
% S-stands for space; T - topology; G - gene; 
mtype = {'ST', 'SG', 'TG', 'S', 'G'}; 
paramMatrix = [1 1 0; 1 0 1; 0 1 1; 1 0 0; 0 0 1]; 

% load empirical data
load('HCPparc20dens_4modelling.mat')
% load genetic data to get the INDs for regions with no data
load('100DS180scaledRobustSigmoidNSGRNAseqQC1Lcortex_ROI_NOdistCorrSurface.mat', 'averageCoexpression'); 
INDnan = find(all(isnan(averageCoexpression),2)); 
INDdata = setdiff(1:size(averageCoexpression,1), INDnan); 

% put data into energy and correlation cells
ENERGY = cell(1,length(mtype)); 
SPTLCORR = cell(1,length(mtype)); 
BestparamNET = cell(1,length(mtype));
allNET = cell(1,length(mtype));

for m=1:length(mtype)
    
    ENERGY{m} = DATA{m}.E;
    SPTLCORR{m} = DATA{m}.Cs;
    BestparamNET{m} = DATA{m}.Bbest; 
    for p=1:length(DATA{m}.E)
        N = zeros(length(INDdata));
        N(DATA{m}.B{p}) = 1;
        N = N+N';
        allNET{m}{p} = N; 
    end

end


for j = 1:length(individualCOORDs)
    DD(:,:,j) = pdist2(individualCOORDs{j},individualCOORDs{j});
end
Dist = mean(DD,3);

% get empirical network
if strcmp(hemi, 'left')
    E = logical(GrFA(1:size(GrFA)/2, 1:size(GrFA)/2));
    D = Dist(1:size(GrFA)/2, 1:size(GrFA)/2); 
    
    % remove non-existing rows and columns from empirical data
    D(INDnan,:) = []; 
    D(:, INDnan) = []; 
    
    E(INDnan,:) = []; 
    E(:, INDnan) = []; 

else
    E = logical(GrFA);
    D = Dist; 
end

degEMP = degrees_und(E);
BestParamNetworksCORR = cell(1,1);
numMOD = length(mtype);
numNET = length(DATA{1}.Ebest);

switch optimiseWhat
    case 'degCorr' % for degree correlations look for max values
        
        for m=1:numMOD
            for n=1:numNET
                
                % incert NaNs where gene data were missing
                networkSEL = DATA{m}.Bbest{n};

                degnetworkSEL = degrees_und(networkSEL);
                degnetworkSEL(INDnan) = NaN; % for regions with missing data give NaN degree
                BestParamNetworksCORR{m}(n) = corr(degEMP', degnetworkSEL', 'rows', 'complete', 'type', 'Spearman');
            end
        end
        
        % plot top 100 highest correlation values
        [S, INDselected] = plotMODviolin(SPTLCORR, 100, mtype, 'highest');
        set(gca,'FontSize',18)
        ylabel('Degree correlation')
        figureName = sprintf('makeFigures/MODELfit_%s_%s_Highestcorrelation_optimise_%s_genes.png',parcellation, hemi, optimiseWhat);
        
        
        % get CORRELATION values for those selected networks
        CS = cell(1,length(mtype));
        for p=1:length(mtype)
            CS{p} = SPTLCORR{p}(INDselected{p});
        end
        
        
    case 'energy' % for KS look for minimal values
        
        [S, INDselected] = plotMODviolin(ENERGY, 100, mtype, 'lowest');
        set(gca,'FontSize',18)
        ylabel('Model fit (KS)')
        figureName = sprintf('makeFigures/MODELfit_%s_%s_Lowestenergy_optimise_%s_genes.png',parcellation, hemi, optimiseWhat);

        % get the correlations based on lowest energy
        CS = cell(1,length(mtype));
        for p=1:length(mtype)
            CS{p} = SPTLCORR{p}(INDselected{p});
        end
        
end


% select best model based on 10000 runs
F = fields(S);
Fbest = F{1}; % fields are ordered based on mean energy
[~,V] = intersect(mtype, Fbest); 
MenergyIND = INDselected{V}(1); 


% select best parameters based on MenergyIND values 
bestPARAM = DATA{V}.P(MenergyIND,:); 
bestPARAM = bestPARAM.*paramMatrix(V,:); 

% plot CDFs
% figure('color','w');
% set(gcf, 'Position', [500 500 500 750])
% CDFs on best parameters
BESTmodel_networks = allNET{V}(INDselected{V}); 
%BESTmodel_networks = BestparamNET{V};
plot_modellingCDF(E, BESTmodel_networks, D, 100);
% save the figure
print(gcf,figureName,'-dpng','-r300');

% best model over 10000 runs
BESTmodel = allNET{V}{MenergyIND}; 

fprintf ('model %s has the "best" network\n', mtype{V})
fprintf ('parameters are eta=%d, gamma=%d, lambda=%d\n', bestPARAM(1), bestPARAM(2), bestPARAM(3))

% get max KS value for main text
mName = erase(mtype{V}, '-');
maxKS = max(S.(mName));
fprintf ('Max KS for best model in box plot is %d\n', maxKS)

degMOD = degrees_und(BESTmodel);
%degMOD(INDnan) = NaN; % replace degree for non-existing regions with NaN; 

[rS,pS] = corr(degMOD', degEMP', 'rows', 'complete', 'type', 'Spearman');
[rP,pP] = corr(degMOD', degEMP', 'rows', 'complete');
CorrLE_ALL = cat(1, CS{:});

whatLine = 'fit';
figure('color','w');
set(gcf, 'Position', [500 500 500 750])

subplot(9,1,1:3);
scatter(degEMP, degMOD, 150,'MarkerEdgeColor',[90,174,97]/255,'MarkerFaceColor',[1 1 1], 'LineWidth',3);
set(gcf, 'renderer', 'painters')
hold on;
switch optimiseWhat
    case 'energy'
        xlim([0 200]);
        ylim([0 200]);
        % plot degree on brain for empirical and mode network
        d = 25;
        tsEMP = [145-d,125-d,105-d];
        tsMOD = [145-d,125-d,105-d];%[80,75,70];
    case 'degCorr'
        xlim([0 150]);
        ylim([0 150]);
        % plot degree on brain for empirical and mode network
        demp = 25;
        % use different thresholds for degre-optimised modelling - the distribution is very narrow
        tsEMP = [145-demp,125-demp,105-demp];
        tsMOD = [65,60,55];
end

switch whatLine
    case 'identity'
        X = linspace(1,200,20);
        plot(X,X,':','LineWidth',4,'Color', [0.2 0.2 0.2]);
    case 'fit'
        h2 = lsline;
        h2.LineWidth = 4;
        h2.LineStyle = ':';
        h2.Color = [0.2, 0.2, 0.2];
end

axis square
xlabel({'Node degree', 'empirical data'});
ylabel({'Node degree', 'best-fitting model'})

set(gca,'fontsize', 20);

% make separate fistograms for G and other models
has_genes = contains(mtype, 'G'); 
%Corr_G = vertcat(CS{Gmod});
%Corr_notG = vertcat(CS{Gmod==0});
for gg=1:length(mtype)
    
    if has_genes(gg)
        bar_col = [90,174,97]/255;
    else
        bar_col = [69,117,180]/255;
    end
    
    subplot(9,1,gg+4);
    histogram(CS{gg}, 20, 'EdgeColor',bar_col,'FaceColor',[1 1 1], 'LineWidth',3);
    xlim([-0.35 0.35]); ylim([0 15]); box off
    set(gca,'fontsize', 14);
    
    if gg~=5
        set(gca,'xticklabel',[])
    end
    
    if gg==3
        ylabel('Frequency', 'FontSize', 20)
    elseif gg==5
        xlabel('Spearman correlation, \rho', 'FontSize', 20)
    end
    
    set(gca,'ytick',[0 15])
    set(gca,'xtick',[-.2 0 .2])
    
    hold on;
end

figureName = sprintf('makeFigures/MODdegree_EMPvsMOD_%s_%s_%s_genes.png', optimiseWhat,parcellation, hemi);
print(gcf,figureName,'-dpng','-r300');

degEMP_plot = nan(180,1); 
degEMP_plot(INDdata) = degEMP; 

plot_hubsSurface_SO('HCP',degEMP_plot,tsEMP, 'inside', 'lh');
figureName = sprintf('makeFigures/hubsSurface_EMP_%s_%s_%s_%s_genes.png', optimiseWhat, parcellation, 'inside', 'lh');
print(gcf,figureName,'-dpng','-r300');

plot_hubsSurface_SO('HCP',degEMP_plot,tsEMP, 'outside', 'lh');
figureName = sprintf('makeFigures/hubsSurface_EMP_%s_%s_%s_%s_genes.png', optimiseWhat, parcellation, 'outside', 'lh');
print(gcf,figureName,'-dpng','-r300');


sides = {'inside'; 'outside'};
hems = {'lh'};
for s=1:2
    side = sides{s};
    for h=1:length(hems)
        hem = hems{h};
        ds_all = nan(1, 180);
        ds_all(INDdata) = degMOD;
        
        plot_hubsSurface_SO('HCP',ds_all,tsMOD, side, 'lh');
        figureName = sprintf('makeFigures/hubsSurface_bestMOD_%s_%s_%s_%s_genes.png', optimiseWhat, parcellation, side, hem);
        print(gcf,figureName,'-dpng','-r300');
        
    end
end
end