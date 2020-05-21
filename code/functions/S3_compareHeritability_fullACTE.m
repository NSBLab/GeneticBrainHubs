function [heritMatrix, nodeData, groupAdjlog] = S3_compareHeritability_fullACTE(parcellation,tractography,plotWhat,weight2,densThreshold,cvMeasure, plotOptions, n)
 % indTOP, mask
if nargin < 8
    n = 500;
end
whatDistribution = plotOptions.whatDistribution; 
colorOut = plotOptions.colorOut; 
colorIn = plotOptions.colIn; 


if strcmp(parcellation, 'HCP') || strcmp(parcellation, 'cust100')
    numSubc = 10; conFile = 'HCPMMP1ANDfslatlas20';
    cortex = [1:180,191:370];
% elseif strcmp(parcellation, 'aparcaseg')
%     numSubc = 7;
%     cortex = [1:34,42:75];
elseif strcmp(parcellation, 'cust250')
    numSubc = 15; conFile = 'custom500ANDfslatlas20';
    cortex = [1:250,266:515];
% elseif strcmp(parcellation, 'Schaefer200')
%     numSubc = 0; conFile = 'Schaefer200_17net';
%     cortex = [1:200];
% elseif strcmp(parcellation, 'Schaefer400')
%     numSubc = 0; conFile = 'Schaefer400_17net';
%     cortex = [1:400];
% elseif strcmp(parcellation, 'Schaefer900') 
%     numSubc = 0; conFile = 'Schaefer900_17net';
%     cortex = [1:900];
end
    
load(sprintf('twinEdges_%s_%s_%s_%s%d.mat', parcellation, tractography, weight2, cvMeasure, round(densThreshold*100))); 
heritFile = sprintf('heritability_onlyACTEnoOUTLIERS_wSATpValsnew_allEdges_twinEdges_%s_%s_%s_%s%d.mat-1.txt', parcellation, tractography, weight2, cvMeasure, round(densThreshold*100)); 
heritabilityACE = importHeritabilityResult(heritFile); 

% count the number of edges for every model
ACTE = sum(heritabilityACE.heritabilityA ~= 0 & ...
           heritabilityACE.heritabilityC ~= 0 & ...
           heritabilityACE.heritabilityT ~= 0 & ...
           heritabilityACE.heritabilityE ~= 0);    
ACE = sum(heritabilityACE.heritabilityA ~= 0 & ...
          heritabilityACE.heritabilityC ~= 0 & ...
          heritabilityACE.heritabilityT == 0 & ...
          heritabilityACE.heritabilityE ~= 0); 
AE = sum(heritabilityACE.heritabilityA ~= 0 & ...
         heritabilityACE.heritabilityC == 0 & ...
         heritabilityACE.heritabilityT == 0 & ...
         heritabilityACE.heritabilityE ~= 0); 
CE = sum(heritabilityACE.heritabilityA == 0 & ...
         heritabilityACE.heritabilityC ~= 0 & ...
         heritabilityACE.heritabilityT == 0 & ...
         heritabilityACE.heritabilityE ~= 0); 
E = sum(heritabilityACE.heritabilityA == 0 & ...
        heritabilityACE.heritabilityC == 0 & ...
        heritabilityACE.heritabilityT == 0 & ...
        heritabilityACE.heritabilityE ~= 0); 
    
fprintf('ACTE model: %d edges, %d perc\n', ACTE, round(100*ACTE/size(heritabilityACE,1),1)); 
fprintf('ACE model: %d edges, %d perc\n', ACE, round(100*ACE/size(heritabilityACE,1),1)); 
fprintf('AE model: %d edges, %d perc\n', AE, round(100*AE/size(heritabilityACE,1),1)); 
fprintf('CE model: %d edges, %d perc\n', CE, round(100*CE/size(heritabilityACE,1),1)); 
fprintf('E model: %d edges, %d perc\n', E, round(100*E/size(heritabilityACE,1),1)); 

numNodes = size(groupAdjlog,1); 
% reshape heritability vector into the matrix for connected edges get indexes on diagonal
heritMatrix = zeros(numNodes,numNodes);
% mask upper half to get indexes of existing links
C = maskuHalf(groupAdjlog); 
% combine values to a vector for reshaping
% in heritability variable 1st column A, 2nd column C, 3rd column E
switch plotWhat
    case 'Afactor'
heritMatrix(C==1) = heritabilityACE.heritabilityA; % assign heritability values to those edges
    case 'Cfactor'
heritMatrix(C==1) = heritabilityACE.heritabilityC; % assign common environmental values
    case 'Efactor'
heritMatrix(C==1) = heritabilityACE.heritabilityE; % assign unique environmental values
    case 'Tfactor'
heritMatrix(C==1) = heritabilityACE.heritabilityT; % assign twin-specific environmental values
end
% make a full matrix
heritMatrixHalf = maskuHalf(heritMatrix); 
heritMatrix = heritMatrix+heritMatrix';
% try remove links that have 0 heritability
% heritMatrix(heritMatrix==0) = NaN; 
%
nodeData = degrees_und(groupAdjlog); 
% make a curve plot for the whole brain
getMaxVal = RichClubHuman_median(groupAdjlog,heritMatrix, nodeData,'right', whatDistribution, colorOut, colorIn); 
%title('Cortex')
ylabel('Median edge heritability')
set(gcf, 'Position', [500 500 750 550])
set(gca,'fontsize', 18);
ylim([0 max(getMaxVal)+0.5*nanstd(getMaxVal)])

% if strcmp(weight2, 'FA')
%     ylim([0.1 0.71])
% elseif strcmp(weight2, 'standard')
%     ylim([0 0.4])
% end
% get top 500 links for the spatial plot
ind(:,1) = find(~isnan(heritMatrixHalf)); 
ind(:,2) = heritMatrix(~isnan(heritMatrixHalf(:))); 

indS = sortrows(ind,-2); 
% 
 indTOP = indS(1:n,:); 
% BOT = indS(indS(:,2)==0,:); 
% indBOT(:,1) = datasample(BOT(:,1), n); 
% indBOT(:,2) = 0; 
% 
% 
% % % put vallues back to matrix
% mask = zeros(size(heritMatrix(:))); 
% mask(indTOP(:,1)) = 2; 
% mask(indBOT(:,1)) = 1; 
% mask = reshape(mask,size(heritMatrix)); 

end

