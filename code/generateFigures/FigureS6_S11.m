%--------------------------------------------------%
% Figure S6-S9, S11-S12: CGE and MPC
%--------------------------------------------------%
function FigureS6_S11()
whatDWI = 'HCP';
weight = 'standard'; %for GenCog'standard';
parc = 'HCP';
op = selectCONmetrics(parc, weight);
numThr = 4;

plotOptions.colIn = [1 1 1];
plotOptions.colorOut = [82 82 82]/255; %[0,69,41]/255; %[.35 .35 .35]; %[4 90 141]/255; %[.45 .45 .45]; %[5 113 176]/255;
plotOptions.whatDistribution = 'histogram';

[coexpData, A, matrices, coordinates, avWeight] = giveConnExp_HCP(parc,op.tract,op.probe,weight,op.brainPart,op.nRem);
[~, ~, M, ~, avWeightFA] = giveConnExp_HCP(parc,op.tract,op.probe,'FA',op.brainPart, 0);

numLC = size(coexpData.averageCoexpression,1);

GrSC = giveMeRichClub(matrices, coordinates, op.groupConn ,op.densThreshold, false, op.cvMeasure, op.consThr);
GrFA = avWeightFA.*logical(GrSC);
groupAdjlog = logical(GrSC);
nodeData = degrees_und(groupAdjlog);

% Define distance based on the distance between regions on the surface - it's more relevant for CGE then connection distance based on the tract length.
distMatr = coexpData.averageDistance;
% select left hemisphere data for connectivity matrix, degree distribution
groupAdjlog = groupAdjlog(1:numLC, 1:numLC);
nodeData = nodeData(1:numLC);

% use exponential fit to correct for distance effect
% Figure S6
CGEmatrix_uncorrected = corr(coexpData.parcelExpression(:,2:end)');
[CGEmatrix, FitCurve, c] = measure_correctDistance(CGEmatrix_uncorrected, coexpData.averageDistance, 'Correlated gene expression');
CGEmatrix(groupAdjlog==0) = NaN;

% plot CGE for different distance bins as violin plots
% Figure S8
[RvsF, FvsP, dataCell, xThresholds] = plot_distanceViolin(CGEmatrix , distMatr, groupAdjlog, nodeData, op.khub, numThr, 'CGE');
figureName = sprintf('makeFigures/CGEdist_%s_%d.png', parc, round(op.densThreshold*100));
print(gcf,figureName,'-dpng','-r600');

% get p-values comparing rich to peripheral links
for k=1:length(dataCell)
    [~, pRP(k),~,statsRP(k)] = ttest2(dataCell{1,k}{1,1},dataCell{1,k}{3,1}, 'Vartype','unequal', 'Tail', 'right');
    tRP(k) = statsRP.tstat;
end

% Figure S8
% for each of cell-specific gene groups, make R/F/P curves.
[cellGenesUNIQUE,cellGroups] = getCellSpecificGenes('all_gene_data.csv');
nG = length(cellGroups); 

% for each cell-specific list of genes find ones in the expression data
listGenes = getCellExpression(cellGenesUNIQUE,cellGroups,coexpData); 

for g=1:length(cellGroups)
    % extract expression measures for each gene group
    groupExp = listGenes{g,3};
    groupCGE = corr(groupExp');
    % correct CGE for distance effects using the "fitCurve" based on all genes
    groupCGE_DC = groupCGE-FitCurve;
    
    % plot curves
    RichClubHuman(groupAdjlog,groupCGE_DC, nodeData, 'right', plotOptions.whatDistribution, plotOptions.colorOut, plotOptions.colIn);
    ylabel({'Mean spatially-corrected', 'correlated gene expression'})
    set(gcf, 'Position', [500 500 750 550])
    set(gca,'fontsize', 18);
    ylim([-0.05 0.18])
    
    % remove "-" from group name for saving the figure;
    groupName = listGenes{g,1};
    isD = contains(groupName,'-');
    
    if isD
        groupName = erase(groupName,'-');
    end
    
    figureName = sprintf('makeFigures/CGEcurves_%sgenes_%s_%d.png', listGenes{g,1}, parc, round(op.densThreshold*100));
    print(gcf,figureName,'-dpng','-r600');
end


% Figure S7
parcs = {'HCP','random500'};
tract = 'iFOD2';

densThresholds = [0.15, 0.2, 0.25; ... % for HCP parcellation
    0.05 0.1, 0.15]; % for custom 500 parcellation

khubs = [90, 105, 120; ...% for HCP parcellation
    35, 70, 100]; % for custom 500 parcellation


for pa=1:length(parcs)
    parc = parcs{pa};
    
    [coexpData, ~, matrices, coordinates] = giveConnExp_HCP(parc,op.tract,op.probe,weight,op.brainPart,op.nRem);
  
    
    for de=1:size(densThresholds,2)
        densThreshold = densThresholds(pa,de);
        khub = khubs(pa,de);
        numLC = size(coexpData.averageCoexpression,1);
        
        GrSC = giveMeRichClub(matrices, coordinates, op.groupConn ,densThreshold, false, op.cvMeasure, op.consThr);
        groupAdjlog = logical(GrSC);
        nodeData = degrees_und(groupAdjlog);
        
        
        % load data containing connection distance
        distMatr = giveConnDistance(parc, tract, groupAdjlog);
        distMatr = maskuHalf(distMatr);
        % select left hemisphere data for connectivity matrix, degree distribution and distance matrix
        distMatr = distMatr(1:numLC, 1:numLC);
        groupAdjlog = groupAdjlog(1:numLC, 1:numLC);
        nodeData = nodeData(1:numLC);
        
        % use exponential fit to correct for distance effect
        CGEmatrix_uncorrected = corr(coexpData.parcelExpression(:,2:end)');
        [CGEmatrix, FitCurve] = measure_correctDistance(CGEmatrix_uncorrected, coexpData.averageDistance, 'CGE', 'exp', false);
        CGEmatrix(groupAdjlog==0) = NaN;
        
        % plot R/F/P lines for CGE
        whatTail = 'right';
        RichClubHuman(groupAdjlog,CGEmatrix, nodeData, whatTail, plotOptions.whatDistribution, plotOptions.colorOut, plotOptions.colIn);
        ylabel({'Mean spatially-corrected', 'correlated gene expression'})
        set(gcf, 'Position', [500 500 750 550])
        set(gca,'fontsize', 18);
        ylim([-0.05 0.18])
        
        figureName = sprintf('makeFigures/CGEcurves_%s_%d.png', parc, densThreshold);
        print(gcf,figureName,'-dpng','-r600');
    end
end


% Figure S9
% make RFP plot using Monash dataset connectome
[coexpData, A, matrices, coordinates, avWeight] = giveConnExp_GenCog('HCP',op.tract,op.probe,op.conW,op.brainPart,op.nRem);
GrSC_GC = giveMeRichClub(matrices, coordinates, op.groupConn, op.densThreshold, false, op.cvMeasure, op.consThr);

CGEmatrix_uncorrected = corr(coexpData.parcelExpression(:,2:end)');
numLC = size(coexpData.averageCoexpression,1);

groupAdjlogGC = logical(GrSC_GC);
nodeDataGC = degrees_und(groupAdjlogGC);

groupAdjlogGC = groupAdjlogGC(1:numLC, 1:numLC);
nodeDataGC = nodeDataGC(1:numLC);

[CGEmatrixGC, FitCurve] = measure_correctDistance(CGEmatrix_uncorrected, coexpData.averageDistance, 'CGE', 'exp', false);
CGEmatrixGC(groupAdjlogGC==0) = NaN;

whatTail = 'right';
RichClubHuman(groupAdjlogGC,CGEmatrixGC, nodeDataGC, whatTail, plotOptions.whatDistribution, plotOptions.colorOut, plotOptions.colIn);
ylabel({'Mean spatially-corrected', 'correlated gene expression'})
set(gcf, 'Position', [500 500 750 550])
set(gca,'fontsize', 18);
ylim([-0.05 0.1])

figureName = sprintf('makeFigures/CGEcurves_%s_%d_%s.png', parc, densThreshold, whatDWI);
print(gcf,figureName,'-dpng','-r600');
end



