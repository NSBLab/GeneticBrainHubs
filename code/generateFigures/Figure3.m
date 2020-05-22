
%--------------------------------------------------% 
% Figure 3
%--------------------------------------------------%
function Figure3()

whatDWI = 'HCP';
weight = 'standard'; %for GenCog'standard'; 
parc = 'HCP';
op = selectCONmetrics(parc, weight); 
numThr = 4;

plotOptions.colIn = [1 1 1];
plotOptions.colorOut = [82 82 82]/255; %[0,69,41]/255; %[.35 .35 .35]; %[4 90 141]/255; %[.45 .45 .45]; %[5 113 176]/255; 
plotOptions.whatDistribution = 'histogram'; 

if strcmp(whatDWI, 'HCP')
    
    [coexpData, A, matrices, coordinates, avWeight] = giveConnExp_HCP(parc,op.tract,op.probe,weight,op.brainPart,op.nRem);
    [~, ~, M, ~, avWeightFA] = giveConnExp_HCP(parc,op.tract,op.probe,'FA',op.brainPart, 0);
    
elseif strcmp(whatDWI, 'GenCog')

    [coexpData, A, matrices, coordinates, avWeight] = giveConnExp_GenCog(parc,op.tract,op.probe,weight,op.brainPart,op.nRem);
    [~, ~, M, ~, avWeightFA] = giveConnExp_GenCog(parc,op.tract,op.probe,'FA',op.brainPart, 0);
    
end

numLC = size(coexpData.averageCoexpression,1); 

[GrSC, mDIST] = giveMeRichClub(matrices, coordinates, op.groupConn ,op.densThreshold, false, op.cvMeasure, op.consThr);
GrFA = avWeightFA.*logical(GrSC); 
groupAdjlog = logical(GrSC); 
nodeData = degrees_und(groupAdjlog);

% Define distance based on the distance between regions on the surface - it's more relevant for CGE then connection distance based on the tract length. 
distMatr = coexpData.averageDistance; 
% select left hemisphere data for connectivity matrix, degree distribution
groupAdjlog = groupAdjlog(1:numLC, 1:numLC); 
nodeData = nodeData(1:numLC); 

% use exponential fit to correct for distance effect
CGEmatrix_uncorrected = corr(coexpData.parcelExpression(:,2:end)');
[CGEmatrix, FitCurve, c] = measure_correctDistance(CGEmatrix_uncorrected, coexpData.averageDistance, 'Correlated gene expression', 'exp', false);
CGEmatrix(groupAdjlog==0) = NaN; 

% plot R/F/P lines for CGE
whatTail = 'right'; 
RichClubHuman(groupAdjlog,CGEmatrix, nodeData, whatTail, plotOptions.whatDistribution, plotOptions.colorOut, plotOptions.colIn); 
ylabel({'Mean spatially-corrected', 'correlated gene expression'})
set(gcf, 'Position', [500 500 750 550])
set(gca,'fontsize', 18);
ylim([-0.05 0.18])

figureName = sprintf('makeFigures/CGEcurves_%s_%d.png', parc,round(op.densThreshold*100));
print(gcf,figureName,'-dpng','-r600');

% plot CGE for each functional module
load('data/modules/cortex_parcel_network_assignments.mat')

[~,netNamesAll,allData,measureModRICH] = compareMeasure_modules(CGEmatrix, netassignments, groupAdjlog, nodeData, op.khub);
resINTERINTRA = plot_moduleViolin(allData, measureModRICH, netNamesAll, netassignments, CGEmatrix, nodeData, groupAdjlog, 'CGE', op.khub); 

figureName = sprintf('makeFigures/CGEmodules_%s_%d.png', parc, round(op.densThreshold*100));
print(gcf,figureName,'-dpng','-r600');

% get cell-specific gene lists
[cellGenesUNIQUE,cellGroups] = getCellSpecificGenes('all_gene_data.csv');
nG = length(cellGroups); 
% add disorder-specific genes to cellGenesUNIQUE and cellGroups variables
lists = {'iq', 'scz'}; 
N = getGWASgenes(coexpData, lists); 
for gg=1:length(N)
    cellGenesUNIQUE{nG+gg} = string(N{gg}'); 
    cellGroups(nG+gg) = lists{gg}; 
end

% for each cell-specific list of genes find ones in the expression data
listGenes = getCellExpression(cellGenesUNIQUE,cellGroups,coexpData); 
[T, geneList, fig] = GCCttest(coexpData, FitCurve, groupAdjlog, nodeData, op.khub, listGenes); 

figureName = sprintf('makeFigures/GCCdistributions_%s_%d.png',  parc, round(op.densThreshold*100));
print(gcf,figureName,'-dpng','-r600');

% Load BigBrain data
yregress = true; 
[BBmean,microDist,distNorm, BBmpc, mp, BBskew] = getBBdata(parc, yregress);

whatTail = 'right';
nodeDeg = degrees_und(GrSC); 
RichClubHuman(logical(GrSC),BBmpc,nodeDeg, whatTail, plotOptions.whatDistribution, plotOptions.colorOut, plotOptions.colIn);
ylabel('Mean microstructural covariance')
set(gcf, 'Position', [500 500 750 550])
set(gca,'fontsize', 18);
ylim([0.6 1.1])

figureName = sprintf('makeFigures/MPCmean_%s_%d.png', parc, round(op.densThreshold*100));
print(gcf,figureName,'-dpng','-r600');
end


