% make lists of genes for enrichment
whatSets = {'POS'; 'NEG'}; 

for i=1:length(whatSets)
    
    % HCP
    [~, genes_ENTREZ_HCP] = extract_eQTLs('HCP', whatSets{i}); 
    genes_ENTREZ_HCP = table(string(genes_ENTREZ_HCP)); 
    
    fileName = sprintf('data/reeqtls/%ssnps/HCP_%s_listGenes_onlyEntrezID.txt', whatSets{i}, whatSets{i}); 
    writetable(genes_ENTREZ_HCP, fileName, 'WriteVariableNames', false)
    
    
    
    % GenCog
    [~, genes_ENTREZ_GenCog] = extract_eQTLs('GenCog', whatSets{i}); 
    genes_ENTREZ_GenCog = table(string(genes_ENTREZ_GenCog)); 
    
    fileName = sprintf('data/reeqtls/%ssnps/GenCog_%s_listGenes_onlyEntrezID.txt', whatSets{i}, whatSets{i}); 
    writetable(genes_ENTREZ_GenCog, fileName, 'WriteVariableNames', false)
    
end

