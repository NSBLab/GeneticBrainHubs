%--------------------------------------------------% 
% Figure 4
%--------------------------------------------------%
clear all;
close all;

hemi = 'left'; %if using degCorr, select only 'left'
parcellation = 'HCP';
optimiseWhat = 'energy';

switch optimiseWhat
    case 'energy'
        load('GenerativeModelling_EnergyOptim.mat')
        
    case 'degCorr'
        load('GenerativeModelling_CorrOptim.mat')
        
end

rowIND = find(contains(Parcname,parcellation));
load('HCPparc20dens_4modelling.mat')

% get empirical network
if strcmp(hemi, 'left')
    E = logical(GrFA(1:size(GrFA)/2, 1:size(GrFA)/2));
    colIND = 1;
else
    E = logical(GrFA);
    colIND = 2;
end

degEMP = degrees_und(E);
BestParamNetworksCORR = cell(3,2);
numMOD = length(mtype)-1;
numNET = length(BestParamNetworksENERGY{1, 1}{1,1});

switch optimiseWhat
    case 'degCorr' % for degree correlations look for max values

        for m=1:numMOD
            for n=1:numNET
                networkSEL = BestParamNetworks{rowIND,colIND}{1,m}{n};
                degnetworkSEL = degrees_und(networkSEL);
                BestParamNetworksCORR{rowIND, colIND}{1,m}(n) = corr(degEMP', degnetworkSEL', 'type', 'Spearman');
            end
        end

        % plot top 100 highest correlation values
        [S, INDselected] = plotMODviolin(SPTLCORR{rowIND,colIND}, 100, mtype, 'highest');
        set(gca,'FontSize',18)
        ylabel('Degree correlation')
        figureName = sprintf('makeFigures/MODELfit_%s_%s_Highestcorrelation_optimise_%s.png',parcellation, hemi, optimiseWhat);
        print(gcf,figureName,'-dpng','-r600');

        % get CORRELATION values for those selected networks
        CS = cell(1,length(mtype)-1);
        for p=1:length(mtype)-1
            CS{p} = SPTLCORR{rowIND,colIND}{p}(INDselected{p});
        end

        
    case 'energy' % for KS look for minimal values

        [S, INDselected] = plotMODviolin(ENERGY{rowIND,colIND}, 100, mtype, 'lowest');
        set(gca,'FontSize',18)
        ylabel('Model fit (KS)')
        figureName = sprintf('makeFigures/MODELfit_%s_%s_Lowestenergy_optimise_%s.png',parcellation, hemi, optimiseWhat);
        print(gcf,figureName,'-dpng','-r600');
        
        
        % get the correlations based on lowest energy
        CS = cell(1,length(mtype)-1);
        for p=1:length(mtype)-1
            CS{p} = SPTLCORR{rowIND,colIND}{p}(INDselected{p});
        end
        
end


% select best model based on 10000 runs
for t=1:length(mtype)-1
    switch optimiseWhat
        case 'energy'
    Menergy(t) = min(ENERGY{rowIND,colIND}{1,t}(:));
    [~,V] = min(Menergy);
        case 'degCorr'
    Menergy(t) = max(SPTLCORR{rowIND,colIND}{1,t}(:));
    [~,V] = max(Menergy);
    end     
end

BESTmodel = BestNetwork{rowIND,colIND}{1,V};
fprintf ('model %s has the "best" network\n', mtype{V})

% get max KS value for main text
mName = erase(mtype{V}, '-');
maxKS = max(S.(mName));
fprintf ('Max KS for best model in box plot is %d\n', maxKS)

degMOD = degrees_und(BESTmodel);

[rS,pS] = corr(degMOD', degEMP', 'type', 'Spearman');
[rP,pP] = corr(degMOD', degEMP');
CorrLE_ALL = cat(1, CS{:});

whatLine = 'fit'; 
figure('color','w');
set(gcf, 'Position', [500 500 500 750])

subplot(2,1,1);
scatter(degEMP, degMOD, 150,'MarkerEdgeColor',[69,117,180]/255,'MarkerFaceColor',[1 1 1], 'LineWidth',3);
hold on;
switch optimiseWhat
    case 'energy'
xlim([0 200]);
ylim([0 200]);
% plot degree on brain for empirical and mode network
d = 25;
tsEMP = [145-d,125-d,105-d];
tsMOD = [145-d,125-d,105-d];
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

subplot(2,1,2);
histogram(CorrLE_ALL, 20, 'EdgeColor',[69,117,180]/255,'FaceColor',[1 1 1], 'LineWidth',3);
axis square
ylabel('Frequency')
xlabel('Spearman correlation, \rho')
xlim([-0.35 0.35]); box off
set(gca,'fontsize', 20);

figureName = sprintf('makeFigures/MODdegree_EMPvsMOD_%s_%s_%s.png', optimiseWhat,parcellation, hemi);
print(gcf,figureName,'-dpng','-r600');


plot_hubGroupsSurface('HCP',degEMP(1:180),tsEMP, 'inside', 'lh');
figureName = sprintf('makeFigures/hubsSurface_EMP_%s_%s_%s_%s.png', optimiseWhat, parcellation, 'inside', 'lh');
print(gcf,figureName,'-dpng','-r1200');

plot_hubGroupsSurface('HCP',degEMP(1:180),tsEMP, 'outside', 'lh');
figureName = sprintf('makeFigures/hubsSurface_EMP_%s_%s_%s_%s.png', optimiseWhat, parcellation, 'outside', 'lh');
print(gcf,figureName,'-dpng','-r1200');


sides = {'inside'; 'outside'};
hems = {'lh'};
for s=1:2
    side = sides{s};
    for h=1:length(hems)
        hem = hems{h};
        if strcmp(hem, 'lh')
            ds = degMOD(1:180);
        elseif strcmp(hem, 'rh')
            ds = degMOD(181:360);
        end
        
        plot_hubGroupsSurface('HCP',ds,tsMOD, side, 'lh');
        figureName = sprintf('makeFigures/hubsSurface_bestMOD_%s_%s_%s_%s.png', optimiseWhat, parcellation, side, hem);
        print(gcf,figureName,'-dpng','-r1200');
        
    end
end
