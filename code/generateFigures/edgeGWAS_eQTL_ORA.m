% perform over-representation analysis for eQTL genes
% load a list of genes from eQTL mapping
function [pORA_eQTL_HCP, pORA_eQTL_Monash,pORA_eQTL_Consensus, uGenesHCP, uGenesMonash, uGenesConsensus] = edgeGWAS_eQTL_ORA()

uGenesHCP = importGeneList('HCP_listGenes_onlyEntrezID.txt');
uGenesHCP = uGenesHCP(~isnan(uGenesHCP));
THCP = table(uGenesHCP); 
% save consensus genes as a .txt file
writetable(THCP, 'data/reeqtls/HCP_hubeQTL.txt', 'WriteVariableNames',false); 

uGenesMonash = importGeneList('GenCog_listGenes_onlyEntrezID.txt');
uGenesMonash = uGenesMonash(~isnan(uGenesMonash));

TMonash = table(uGenesMonash); 
% save consensus genes as a .txt file
writetable(TMonash, 'data/reeqtls/Monash_hubeQTL.txt', 'WriteVariableNames',false); 

% get gene names for overlapping genes
HCPgeneID = importGENEIDfile('HCP_listGenes_entrezID.csv');
uGenesConsensus = intersect(uGenesHCP, uGenesMonash);
Tconsensus = table(uGenesConsensus); 
% save consensus genes as a .txt file
writetable(Tconsensus, 'data/reeqtls/consensusHCPMonash.txt', 'WriteVariableNames',false); 

[~, indSEL] = intersect(HCPgeneID.converted_alias, uGenesConsensus);
gNAMES = HCPgeneID.name(indSEL);

% compare the overlap between lists
% select only annotated genes

whatAnnotation = 'PSYCHENCODE'; 

n = length(uGenesHCP); % genes in group 1
x = length(intersect(uGenesHCP, uGenesMonash)); % genes between groups
D = length(uGenesMonash); % genes in group 2
switch whatAnnotation
    case 'GTEx'
        N = 15626; % total number of proteing-coding genes considered in eQTL mapping.
    case 'PSYCHENCODE'
        N = 25699; % total number of genes considered in eQTL mapping.
end

% calculate the complement of the hypergeometric pdf
% This will give a probability of getting more than x gene overlap
p = hygecdf(x,N,D,n, 'upper');
fprintf('The probability of overlap more than %d, is p=%d\n', x, p)

% compare overlap with disorder and IQ GWAS lists of genes
% choose eQTL annotation used: PSYCHENCODE or GTEx; 

pORA_eQTL_Monash = eQTL_ORA(uGenesMonash, N, whatAnnotation);

pORA_eQTL_HCP = eQTL_ORA(uGenesHCP, N, whatAnnotation); 

pORA_eQTL_Consensus = eQTL_ORA(uGenesConsensus, N, whatAnnotation); 

end


