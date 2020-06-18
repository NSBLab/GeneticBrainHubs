% ---------------------------------------------------------------------
% Find if genes that are outliers in non-cell specific list are
% ---------------------------------------------------------------------
function pORA = eQTL_ORA(eQTLgenes, M, whatAnnotation)
% test if disease-related genes are over-represented in the list of eQTL
% mapped genes
lists = {'scz', 'adhd', 'aut', 'bip', 'mdd', 'iq'};
% options specific for GWAS-based lists
brainParts = {'Cortex'}; %, 'Frontal_Cortex_BA9', 'Anterior_cingulate_cortex_BA24'};
numTissues = length(brainParts);

for l=1:length(lists)
    whatGeneSet = lists{l};
    
    % when HAR genes are selected
    if strcmp(whatGeneSet, 'har')
        % list of HAR genes was meade using getHARgenes.m, these were converted to entrezIDs
        % link to gene name-to ID conversion here: https://biit.cs.ut.ee/gplink/l/1t5sZyTNRF
        listHAR = importHARgeneIDs('gProfiler_hsapiens_18-06-2020_11-22-26.csv');
        listHAR(isnan(listHAR.converted_alias),:) = [];
        results.(whatGeneSet) = listHAR.converted_alias;
        resultsNames.(whatGeneSet) = listHAR.initial_alias;
    else
        % for non-har genes import based on annotation
        if strcmp(whatAnnotation, 'PSYCHENCODE')
            [results,resultsNames] = selectGWASgenes_eQTL(whatGeneSet, brainParts, numTissues, 'oneList', 1, whatAnnotation);
        else
            [results] = selectGWASgenes_eQTL(whatGeneSet, brainParts, numTissues, 'oneList', 1, whatAnnotation);
        end
    end

    if sum(strcmp(fieldnames(results), whatGeneSet)) == 1
        selectedGenes = results.(whatGeneSet);
        if strcmp(whatAnnotation, 'PSYCHENCODE')
            selectedGenesNames = resultsNames.(whatGeneSet);
        end
    end
        

    
    if ~isempty(results.(whatGeneSet))
        % X - number of genes that need to exceed
        % what is the actual number of genes from a list among outliers
        [empOverlap,~,IND] = intersect(eQTLgenes, selectedGenes, 'stable');
        % X = 0:length(empOverlap);
        % M - total population size (all genes that were considered)
        % K - number of items with a desired characteristic (number of genes in a
        % specific list)
        K = length(selectedGenes);
        % N - number of genes selected (number of genes in eQTL list)
        N = length(eQTLgenes);
        if ~isempty(empOverlap)
        p = hygecdf(length(empOverlap),M,K,N,'upper');
        pORA.(whatGeneSet).p = p;
        pORA.(whatGeneSet).Noverlap = length(empOverlap);
        pORA.(whatGeneSet).Ndisorder = length(selectedGenes);
        pORA.(whatGeneSet).Genes = empOverlap;
        if strcmp(whatAnnotation, 'PSYCHENCODE')
            pORA.(whatGeneSet).GeneNames = selectedGenesNames(IND);
        end
        end
    end
    clearvars selectedGenes
end
end
