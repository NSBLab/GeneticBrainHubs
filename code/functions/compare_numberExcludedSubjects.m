function [RvsF, FvsP, dataCell] = compare_numberExcludedSubjects(plotWhat)

if nargin < 1
    plotWhat = 'BOTH';
end

parcellation = 'HCP';
weight = 'FA';
op = selectCONmetrics(parcellation, weight);
numThr = 4;

% load connectivity data
load('twinEdges_HCP_iFOD2_FA_strength20.mat')
V = readtable('subjectsRemoved_allEdges_twinEdges_HCP_iFOD2_FA_strength20.mat.txt');
H = readtable('heritabilityACTEnoOUTLIERSnew_wSATpVals_allEdges_twinEdges_HCP_iFOD2_FA_strength20.mat-1.txt'); 


nodeData = degrees_und(groupAdjlog);
numNodes = size(groupAdjlog,1);
% reshape heritability vector into the matrix for connected edges get indexes on diagonal
valMatrix = zeros(numNodes, numNodes);

% mask upper half to get indexes of existing links
C = maskuHalf(groupAdjlog);

switch plotWhat
    case 'DZ'
        valMatrix(C==1) = V.heritabilitySubRemDZ;
        selMeasure = V.heritabilitySubRemDZ;
    case 'MZ'
        valMatrix(C==1) = V.heritabilitySubRemMZ;
        selMeasure = V.heritabilitySubRemMZ; 
    case 'BOTH'
        valMatrix(C==1) = V.heritabilitySubRemMZ+V.heritabilitySubRemDZ;
        selMeasure = V.heritabilitySubRemMZ+V.heritabilitySubRemDZ;
    case 'VARrem'
        valMatrix(C==1) = V.edgeVarREM;
        selMeasure = V.edgeVarREM;
    case 'VARorig'
        valMatrix(C==1) = V.edgeVarORIG;
        selMeasure = V.edgeVarORIG;
end

if strcmp(plotWhat, 'BOTH')
    % mean number of sujects excluded
    num_median = median(selMeasure);
    num_max = max(selMeasure);
    
    %total number of subjcts
    num_tot = sum(sum(~isnan(Output_MZ(:,1:3,1))) + sum(~isnan(Output_DZ(:,1:3,1))));
    
    % proc of excluded subjects
    perc_median = num_median*100/num_tot;
    perc_max = num_max*100/num_tot;
end

valMatrix = valMatrix+valMatrix';

plotOptions.colIn = [1 1 1];
plotOptions.colorOut = [82 82 82]/255; %[0,69,41]/255; %[.35 .35 .35]; %[4 90 141]/255; %[.45 .45 .45]; %[5 113 176]/255;
plotOptions.whatDistribution = 'histogram';

heritMatrixHalf = zeros(numNodes, numNodes);
heritMatrixHalf(C==1) = H.heritabilityA; 
heritMatrix = heritMatrixHalf+heritMatrixHalf';
heritMatrixHalf = maskuHalf(heritMatrix);

[RvsF, FvsP, dataCell,xThresholds,f0] = plot_distanceViolin(heritMatrixHalf, valMatrix, groupAdjlog, nodeData, op.khub, numThr, 'Heritability');

end
