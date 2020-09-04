function [pORA_eQTL, genes_ENTREZ] = extract_eQTLs(whatDATA, whatSET)

if nargin <1
    whatDATA = 'HCP'; 
end
if nargin <2
    whatSET = 'POS'; 
end

% 1. SNP-gene mapping based on the annotation file
fileClumped = sprintf('%s_list_clumpedSNPs_p1_01_p2_05_withRS_FINAL.txt', whatDATA); 
SNPIDs = readtable(fileClumped, 'ReadVariableNames',false); 

% 2. load relevant SNP file
targetSNP = readtable(sprintf('%s_SNPs%s_RvsP_beta_10perc_clumped.txt', whatDATA, whatSET)); 

% 3. get their rs IDs based on the list of clumped SNPs
[~, IND] = intersect(SNPIDs.Var2, targetSNP.SNPs);
targetSNPrs = SNPIDs.Var1(IND); 

% 4. find genes based on annotation from PSYCHencode
genes = importpec_genes('pec_genes.annot_FINAL.txt'); 
genesONLY = importpec_genesONLY('pec_genes.annot_FINAL.txt'); 

GIND = cell(size(targetSNPrs,1),1); 
for i=1:size(targetSNPrs)
    
    % add spaces around rsID so only the complete matches are found
    c = strcat({' '},targetSNPrs{i}, {' '}); 
    GIND{i} = find(contains(genes, c{:})); 
    
end

% get all genes together
G = unique(vertcat(GIND{:}));
genesSEL = genesONLY(G); 

% 5. take the overlap with DER-08a_hg19_eQTL.significant.txt - these are significant eQTLs
Gsig = readtable('DER-08a_hg19_eQTL.significant.txt'); 
Gsig_names = unique(extractBefore(Gsig.gene_id,'.')); 
genesSIG = intersect(genesSEL, Gsig_names); 

% 6. find their entrezIDs, only genes with entrezIDs are used [those should have a defined function and also needed in the enrichment]
% read in BIOMART - gene ID matching
BIOMART = readtable('BIOMART_geneIDs.txt'); 
[~, IND_BIO] = intersect(BIOMART.ensembl_gene_id, genesSIG); 
genes_ENTREZ = BIOMART.entrezgene_id(IND_BIO); 
genes_ENTREZ(contains(genes_ENTREZ, 'NA')) = [];
genes_ENTREZ = str2double(genes_ENTREZ); 

% 7. test ORA
whatAnnotation = 'PSYCHENCODE'; 
N = length(genesONLY); % total number of genes considered in eQTL mapping in the original, N = 25699 is hard-coded; 
pORA_eQTL = eQTL_ORA(genes_ENTREZ, N, whatAnnotation); 

end

