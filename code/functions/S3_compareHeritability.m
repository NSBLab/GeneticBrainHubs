function [heritMatrix, nodeData, groupAdjlog, mask, data_export] = S3_compareHeritability(parcellation,tractography,plotWhat,weight2,densThreshold,cvMeasure, plotOptions, onlyACTE,n)
 % indTOP, mask
if nargin < 8
    n = 100;
    onlyACTE=false; 
end

if nargin < 9
    n = 100;
end

whatDistribution = plotOptions.whatDistribution; 
colorOut = plotOptions.colorOut; 
colorIn = plotOptions.colIn; 


if ~onlyACTE

heritFile = sprintf('heritabilityACTEnoOUTLIERSnew_wSATpVals_allEdges_twinEdges_%s_%s_%s_%s%d.mat-1.txt', parcellation, tractography, weight2, cvMeasure, round(densThreshold*100)); 
%heritFile = sprintf('heritabilityACTEwithOUTLIERS_wSAT_allEdges_twinEdges_%s_%s_%s_%s%d.mat.txt', parcellation, tractography, weight2, cvMeasure, round(densThreshold*100)); 
heritabilityACE = importHeritabilityResultwP(heritFile); 

else 
    
heritFile = sprintf('heritability_onlyACTEnoOUTLIERSnew_wSATpVals_allEdges_twinEdges_%s_%s_%s_%s%d.mat-1.txt', parcellation, tractography, weight2, cvMeasure, round(densThreshold*100)); 
heritabilityACE = importHeritabilityResultwP(heritFile); 
end
    
load(sprintf('twinEdges_%s_%s_%s_%s%d.mat', parcellation, tractography, weight2, cvMeasure, round(densThreshold*100))); 


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
    case 'Avariance'
        
end
% make a full matrix
heritMatrixHalf = maskuHalf(heritMatrix); 
heritMatrix = heritMatrix+heritMatrix';
%
nodeData = degrees_und(groupAdjlog); 
% make a curve plot for the whole brain 
[getMaxVal,data_export] = RichClubHuman(groupAdjlog,heritMatrix, nodeData,'right', whatDistribution, colorOut, colorIn);
%getMaxVal = RichClubHuman_median(groupAdjlog,heritMatrix, nodeData,'right', whatDistribution, colorOut, colorIn);
ylabel('Mean edge heritability')
set(gcf, 'Position', [500 500 750 550])
set(gca,'fontsize', 20);
ylim([0 max(getMaxVal)+0.5*nanstd(getMaxVal)])

% get top 500 links for the spatial plot if Afactor is selected
if strcmp(plotWhat, 'Afactor')
    ind(:,1) = find(~isnan(heritMatrixHalf));
    ind(:,2) = heritMatrix(~isnan(heritMatrixHalf(:)));
    
    indS = sortrows(ind,-2);
    indTOP = indS(1:n,:);
    BOT = indS(indS(:,2)==0,:);
    indBOT(:,1) = datasample(BOT(:,1), n);
    indBOT(:,2) = 0;
    
    
    % put vallues back to matrix
    mask = zeros(size(heritMatrix(:)));
    mask(indTOP(:,1)) = 2;
    mask(indBOT(:,1)) = 1;
    mask = reshape(mask,size(heritMatrix));
else
    mask = [];
    
end

end

