% perform over-representation analysis for eQTL genes
% load a list of genes from eQTL mapping
function [pORA_eQTL_Monash, pORA_eQTL_HCP] = edgeGWAS_eQTL_ORA()

uGenesHCP = importGeneList('HCP_listGenes_onlyEntrezID.txt');
uGenesHCP = uGenesHCP(~isnan(uGenesHCP));

uGenesMonash = importGeneList('GenCog_listGenes_onlyEntrezID.txt');
uGenesMonash = uGenesMonash(~isnan(uGenesMonash));

% compare the overlap between lists
% select only annotated genes

n = length(uGenesHCP); % genes in group 1
x = length(intersect(uGenesHCP, uGenesMonash)); % genes between groups
D = length(uGenesMonash); % genes in group 2
N = 25699; % total number of proteing-coding genes considered in eQTL mapping.

% calculate the complement of the hypergeometric pdf
% This will give a probability of getting more than x gene overlap
p = hygecdf(x,N,D,n, 'upper');
fprintf('The probability of overlap more than %d, is p=%d\n', x, p)

% get gene names for overlapping genes
HCPgeneID = importGENEIDfile('HCP_listGenes_entrezID.csv');
gENTREZID = intersect(uGenesHCP, uGenesMonash);
[~, indSEL] = intersect(HCPgeneID.converted_alias, gENTREZID);
gNAMES = HCPgeneID.name(indSEL);

% compare overlap with disorder and IQ GWAS lists of genes
whatAnnotation = 'PSYCHENCODE'; % 'GTEx'; 
pORA_eQTL_Monash = eQTL_ORA(uGenesMonash, N, whatAnnotation);

% compare overlap with disorder and IQ GWAS lists of genes
pORA_eQTL_HCP = eQTL_ORA(uGenesHCP, N, whatAnnotation); 
end


