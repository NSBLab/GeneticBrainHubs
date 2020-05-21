function results = selectGWASgenes(coexpData, disorder, brainParts, numTissues, listGenes, correctMultiple)
% This part of the script
% 1. reads in gene lists based on GWAS data for each disorder;
% 2. corrects p values within each list
% 3.
%   a. makes a single list of significant genes across cortical areas (overlap, geneListINT)
%   b. makes a single slit of significant genes across cortical areas (combination of all, geneListALL)
% 4. gets region x gene matrix for each gene list for each disorder


% get list across all tissues
for i=1:numTissues
    
    filename = sprintf('data/emagma_results/%s_emagma_results/emagma_Brain_%s.csv', disorder, brainParts{i});
    emagmaBrain = importEMAGMA(filename);
    
    if correctMultiple
        pthr = 0.05/size(emagmaBrain,1);
    else
        pthr = 0.05;
    end
    
    geneList.(disorder).(brainParts{i}) = emagmaBrain.GENE(emagmaBrain.P<pthr,:);
    
end

switch listGenes
    case 'combined'
        geneListEXP.(disorder) = unique(vertcat(geneList.(disorder).(brainParts{1}), geneList.(disorder).(brainParts{2}), geneList.(disorder).(brainParts{3})));
    case 'intersect'
        geneListEXP.(disorder) = intersect(geneList.(disorder).(brainParts{1}), intersect(geneList.(disorder).(brainParts{2}), geneList.(disorder).(brainParts{3})));
    case 'oneList'
        geneListEXP.(disorder) = geneList.(disorder).(brainParts{i});
end

% find genes from both lists in gene expression data
[~, indSEL] = intersect(coexpData.probeInformation.EntrezID, geneListEXP.(disorder));

if ~isempty(indSEL)
    % extract expression values from region x gene matrix (add +1 as the first column is region number).
    geneReg.(disorder) = coexpData.parcelExpression(:,indSEL+1);
    
    % visualise ordered gene expression matrices for both lists to see the structure cluster expression matrix for visualisation
    % find nan rows
    %         indNAN = find(isnan(geneReg.(disorders{d})(:,1)));
    %         expMatr = geneReg.(disorders{d});
    %         expMatr(indNAN,:) = [];
    %
    %         % reorder by similarity
    %         ordG = BF_ClusterReorder(expMatr');
    %         ordR = BF_ClusterReorder(expMatr);
    %         figure; imagesc(expMatr(ordR,ordG)); title(sprintf('%s in %s list',disorders{d}, listGenes));
    %         xlabel('genes'); ylabel('regions')
    
    if length(indSEL)<10
        [coeff,score,latent,tsquared,explained,mu] = pca(geneReg.(disorder),'NumComponents',length(indSEL));
    else
        [coeff,score,latent,tsquared,explained,mu] = pca(geneReg.(disorder),'NumComponents',10);
    end
    results.PC = score;
    results.PCexplained = explained;
    results.PCloadings = coeff;
    results.geneName = coexpData.probeInformation.GeneSymbol(indSEL,:);
    results.geneID = coexpData.probeInformation.EntrezID(indSEL,:);
    results.geneExp = geneReg.(disorder);
else
    results = struct; 
end

