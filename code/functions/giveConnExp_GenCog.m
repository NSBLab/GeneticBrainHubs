function [coexpData, A, matrices, coordinates, avWeight, subjects] = giveConnExp_GenCog(parc,tract,probeSelection,weight,brainPart,nRem)
                                                                                   
if nargin<5
    weight = 'standard';
    brainPart = 'Lcortex';
end


if strcmp(parc, 'aparcaseg')
    parcellation = 'aparc+aseg';
    LC = 1:34;
    RC = 42:75;

elseif strcmp(parc, 'HCP')
    parcellation = 'HCPMMP1ANDfslatlas20';
    LC = 1:180;
    RC = 191:370;

elseif strcmp(parc, 'cust100')
    parcellation = 'custom200ANDfslatlas20';
    LC = 1:100;
    RC = 111:210;

elseif strcmp(parc, 'cust250')
    parcellation = 'custom500ANDfslatlas20';
    LC = 1:250;
    RC = 261:510;
end

DS = 100;
CON = load(sprintf('GenCog_%s_%s_NOSIFT_%s.mat', parcellation, tract, weight)); 
ADJS = CON.ADJS(~cellfun('isempty',CON.ADJS));

cort = [LC,RC];
if strcmp(brainPart, 'LRcortex')
    keepNodes = cort;
elseif strcmp(brainPart, 'Lcortex')
    keepNodes = LC;
elseif strcmp(brainPart, 'wholeBrain')
    keepNodes = 1:size(ADJS{1},1);
end

subjectsOK = importSubjects('freesurfer_allSubjects.txt'); 
subjects = intersect(subjectsOK,CON.subs); 

fprintf('issues with freesurfer: n=%d\n', length(setdiff(CON.subs, subjectsOK))); 
numNodes = length(keepNodes);
numSubjects = length(subjects);
matrices = cell(1,numSubjects);
coordinates = cell(1,numSubjects);
A = zeros(numNodes,numNodes,numSubjects);
den = zeros(numSubjects,1); 
for s=1:numSubjects
    m = find(CON.subs==subjects(s)); 
    % remove 25% of weakest links
    ADJS{m}(ADJS{m}<nRem) = 0; % remove less than 10 streamlines
    if ~strcmp(brainPart, 'wholeBrain')
        
        matrices{s} = ADJS{m}(keepNodes, keepNodes);
        A(:,:,s) = matrices{s};
        coordinates{s} = CON.coords{m}(keepNodes,:); % assigne the same coords for all subjects

    else
        
        matrices{s} = ADJS{m};
        A(:,:,s) = matrices{s};
        coordinates{s} = CON.coords{m}; % assigne the same coords for all subjects


    end
    den(s) = density_und(matrices{s}); 
    
end

subRem = find(den <= mean(den)-3*std(den)); 
fprintf('issues with density: n=%d\n', length(subRem)); 

A(:,:,subRem) = []; 
matrices(subRem) = []; 
coordinates(subRem) = []; 
subjects(subRem) = []; 

nanA = A;
nanA(nanA==0) = NaN;
avWeight = nanmean(nanA,3);
avWeight(isnan(avWeight)) = 0;
% load expression data
coexpData = load(sprintf('%dDS%dscaledRobustSigmoidNSG%sQC1Lcortex_ROI_NOdistCorrSurface.mat', DS, max(LC), probeSelection));
end