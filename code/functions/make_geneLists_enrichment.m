% make lists of genes for enrichment
whatSets = {'POS'; 'NEG'}; 

for i=1:length(whatSets)
    
    [~, genes_ENTREZ_HCP] = extract_eQTLs('HCP', whatSets{i}); 
    fileName = sprintf('data/reeqtls/%ssnps/HCP_%s_listGenes_onlyEntrezID.txt', whatSets{i}, whatSets{i}); 
    dlmwrite(fileName,genes_ENTREZ_HCP)
    
    [~, genes_ENTREZ_GenCog] = extract_eQTLs('GenCog', whatSets{i}); 
    fileName = sprintf('data/reeqtls/%ssnps/GenCog_%s_listGenes_onlyEntrezID.txt', whatSets{i}, whatSets{i}); 
    dlmwrite(fileName,genes_ENTREZ_GenCog)
    
end

