function [coexpData, A, matrices, coordinates, avWeight, SUBS] = giveConnExp_HCP(parc,tract,probeSelection, weight,brainPart, nRem)

if nargin<5
    weight = 'standard';
    brainPart = 'Lcortex';
end


if strcmp(parc, 'aparcaseg')
    parcellation = 'aparc+aseg';
    LC = 1:34;
    RC = 42:75;
    nROIs = 82;
elseif strcmp(parc, 'HCP')
    parcellation = 'HCPMMP1ANDfslatlas20';
    LC = 1:180;
    RC = 191:370;
    nROIs = 360;
elseif strcmp(parc, 'cust100') || strcmp(parc, 'cust200') || strcmp(parc, 'custom200') ||  strcmp(parc, 'custom100')
    parcellation = 'custom200ANDfslatlas20';
    LC = 1:100;
    RC = 111:210;
    nROIs = 220;
elseif strcmp(parc, 'cust250') || strcmp(parc, 'cust500') || strcmp(parc, 'custom500') ||  strcmp(parc, 'custom250')
    parcellation = 'custom500ANDfslatlas20';
    LC = 1:250;
    RC = 261:510;
    nROIs = 520;
elseif strcmp(parc, 'random200')
    parcellation = 'random200ANDfslatlas20';
    LC = 1:100;
    RC = 111:210;
    nROIs = 220;
elseif strcmp(parc, 'random500')
    parcellation = 'random500ANDfslatlas20';
    LC = 1:250;
    RC = 261:510;
    nROIs = 520;
elseif strcmp(parc, 'Schaefer400')
    parcellation = 'Schaefer400_17net_acpc';
    LC = 1:200;
    RC = 201:400;
    nROIs = 400; 
end


load(sprintf('HCP_%s_%s_NOSIFT_%s.mat', parcellation, tract, weight))
        % remove subject 299
ADJS{299} = []; ADJS = ADJS(~cellfun('isempty',ADJS));
COG{299} = []; COG = COG(~cellfun('isempty',COG));
SUBS(299) = [];


%identify subjects that have density 2DS lower than average
% for su=1:length(ADJS)
%     ds(su) = density_und(ADJS{su});
% end
% indREM = find(ds< (mean(ds)-2*std(ds)) ); 
% ADJS(indREM) = []; 
% COG(indREM) = []; 
% SUBS(indREM) = []; 

DS = 100;

% stack matrices to 3D matrix
% if strcmp(parcellation, 'aparc+aseg')
%     LC = 1:34;
%     RC = 42:75;
%     nROIs = 82;
% elseif strcmp(parcellation, 'HCPMMP1ANDfslatlas20')
%     LC = 1:180;
%     RC = 191:370;
%     nROIs = 360;
% elseif strcmp(parcellation, 'custom200ANDfslatlas20')
%     LC = 1:100;
%     RC = 111:210;
%     nROIs = 220;
% end

coexpData = load(sprintf('%dDS%dscaledRobustSigmoidNSG%sQC1Lcortex_ROI_NOdistCorrSurface.mat', DS, max(LC), probeSelection));

cort = [LC,RC];
if strcmp(brainPart, 'LRcortex')
    keepNodes = cort;
elseif strcmp(brainPart, 'Lcortex')
    keepNodes = LC;
elseif strcmp(brainPart, 'wholeBrain')
    keepNodes = 1:size(ADJS{1},1);
end

numNodes = length(keepNodes);
numSubjects = length(ADJS);
matrices = cell(1,size(ADJS,2));
coordinates = cell(1,size(ADJS,2));
A = zeros(numNodes,numNodes,numSubjects);
%pRem = 0.50;    


for m=1:length(ADJS)
    % remove 25% of weakest links
    ADJS{m}(ADJS{m}<nRem) = 0; % remove less than 10 streamlines
    if ~strcmp(brainPart, 'wholeBrain')
        matrices{m} = ADJS{m}(keepNodes, keepNodes);
        A(:,:,m) = matrices{m};
        coordinates{m} = COG{m}(keepNodes,:); % assigne the same coords for all subjects
    else
        matrices{m} = ADJS{m};
        A(:,:,m) = matrices{m};
        coordinates{m} = COG{m};
    end
    
end

nanA = A;
nanA(nanA==0) = NaN;
avWeight = nanmean(nanA,3);
avWeight(isnan(avWeight)) = 0;
end