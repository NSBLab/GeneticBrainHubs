function [geneListEXP,geneListNames] = selectGWASgenes_eQTL(disorder, brainParts, numTissues, listGenes, correctMultiple, whatAnnotation)
% This part of the script
% 1. reads in gene lists based on GWAS data for each disorder;
% 2. corrects p values within each list
% 3. makes a single list of significant genes across cortical areas


% get list across all tissues
for i=1:numTissues
    
    switch whatAnnotation
        case 'GTEx'
    filename = sprintf('data/emagma_results/%s_emagma_results/emagma_Brain_%s.csv', disorder, brainParts{i});
    emagmaBrain = importEMAGMA(filename);
        case 'PSYCHENCODE'
    filename = sprintf('data/emagma_results/emagma_pec/pec_%s.genes.out.annot.csv', disorder);
    emagmaBrain = importEMAGMApec(filename); 
    emagmaBrain.GENE = emagmaBrain.ENTREZID; 
    % remove NaN entries
    emagmaBrain(isnan(emagmaBrain.GENE),:) = []; 
    end
    
    if correctMultiple
        pthr = 0.05/size(emagmaBrain,1);
    else
        pthr = 0.05;
    end
    
    geneList.(disorder).(brainParts{i}) = emagmaBrain.GENE(emagmaBrain.P<pthr,:);
    if strcmp(whatAnnotation, 'PSYCHENCODE')
    geneListNames.(disorder) = emagmaBrain.SYMBOL(emagmaBrain.P<pthr,:);
    end
    
end

switch listGenes
    case 'combined'
        geneListEXP.(disorder) = unique(vertcat(geneList.(disorder).(brainParts{1}), geneList.(disorder).(brainParts{2}), geneList.(disorder).(brainParts{3})));
    case 'intersect'
        geneListEXP.(disorder) = intersect(geneList.(disorder).(brainParts{1}), intersect(geneList.(disorder).(brainParts{2}), geneList.(disorder).(brainParts{3})));
    case 'oneList'
        geneListEXP.(disorder) = geneList.(disorder).(brainParts{i});
 
end

end

