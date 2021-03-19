function [heritMatrix, nodeData, groupAdjlog, realTrajectory_save] = S3_compareHeritability_onlytwin(plotOptions, plotWhat)
 % indTOP, mask
 if nargin <2
     plotWhat = 'Afactor';
 end
 

whatDistribution = plotOptions.whatDistribution; 
colorOut = plotOptions.colorOut; 
colorIn = plotOptions.colIn; 

heritFile = 'heritabilityACEnoOUTLIERSnew_wSATpVals_allEdges_variance_twinEdges_HCP_iFOD2_FA_strength20_only_twin.mat.txt'; 
heritabilityACE = importHeritabilityResultwP(heritFile); 
load('twinEdges_HCP_iFOD2_FA_strength20.mat', 'groupAdjlog')
% count the number of edges for every model    
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
%
nodeData = degrees_und(groupAdjlog); 
% make a curve plot for the whole brain 
RichClubHuman(groupAdjlog,heritMatrix, nodeData,'right', whatDistribution, colorOut, colorIn);
set(gcf, 'Position', [500 500 750 550])
set(gca,'fontsize', 20);


end

