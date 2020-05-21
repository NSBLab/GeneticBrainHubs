% visualise different ranges of heritability values
clear all;

parcellation = 'HCP';
tractography = 'iFOD2';
conWeight = 'FA'; 
densThreshold = 0.2; 
cvMeasure = 'strength'; 
numThr = 4;
khub = 105; 

plotOptions.colIn = [1 1 1];
plotOptions.colorOut = [4 90 141]/255;
plotOptions.whatDistribution = 'histogram'; 

[heritMatrix, nodeData, groupAdjlog, indTOP, mask] = S3_compareHeritability(parcellation,tractography,conWeight,densThreshold,cvMeasure, plotOptions); 
load('cortex_parcel_network_assignments.mat')

netNamesAll = {'VIS1';'VIS2';'SM';'CO';'DAN';'LAN';'FPN';'AUD';'DMN';'PM';'VM';'OA'}; 

T0 = 0; 
TL = 0.1; 
TH = 0.5; 

% edges with 0 heritability
maskH0 = heritMatrix==T0; 
vals{1,1} = groupAdjlog.*maskH0; 
vals{1,2} = sprintf('H=%d, %d edges', T0, sum(vals{1,1}(:)/2));



% edges with very low heritability
maskHL = (heritMatrix>T0 & heritMatrix<TL);
vals{2,1} = groupAdjlog.*maskHL; 
vals{2,2} = sprintf('H<%d, %d edges', TL, sum(vals{2,1}(:)/2)); 

% edges with high heritability
maskHH = heritMatrix>TH;
vals{3,1} = groupAdjlog.*maskHH; 
vals{3,2} = sprintf('H>%d, %d edges', TH, sum(vals{3,1}(:)/2)); 

numVals = 200; 
T = table; 

for i=1:length(vals)
   
    vals{i,3} = find(triu(vals{i,1})); 
    y = datasample(vals{i,3},numVals, 'Replace',false);
    M = zeros(size(groupAdjlog));

    [row,col] = ind2sub(size(groupAdjlog,1),y); 
    for k=1:numVals
            M(row(k),col(k)) = 1;
    end
    
    fileName = sprintf('data/testing/herit_%d_numEdges%d.txt', i, numVals);
    dlmwrite(fileName, M, 'delimiter', '\t')

    figure; 
    set(gcf, 'Position', [10 10 1800 600]);
    
    subplot(1,2,1); 
    imagesc(vals{i,1}); colormap(parula);
    title(vals{i,2}); 
    
    subplot(1,2,2); 
    [out outPC outNorm] = plotClassifiedEdges(vals{i,1}, netassignments, 3, netNamesAll);
    title(vals{i,2}); 
    
end



