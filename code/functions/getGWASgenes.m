% select data for GWAS eQTLs
function N = getGWASgenes(coexpData, lists)

brainParts = {'Cortex'}; 
numTissues = length(brainParts);

for l=1:length(lists)
    whatGeneSet = lists{l};
    results = selectGWASgenes_eQTL(whatGeneSet, brainParts, numTissues, 'oneList', 1);
    [~, ind] = intersect(coexpData.probeInformation.EntrezID, results.(whatGeneSet), 'stable'); 
    % get gene names
    N{l} = coexpData.probeInformation.GeneSymbol(ind); 
end
end