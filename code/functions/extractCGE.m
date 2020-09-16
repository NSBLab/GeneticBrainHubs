% extract CGE data for :
% 1. all genes
% 2. genes from enriched categories;
% 3. save BigBrain matrix - this one has NaNs
function extractCGE()

weight = 'standard';
parc = 'HCP';
op = selectCONmetrics(parc, weight);

[coexpData, A, matrices, coordinates, avWeight] = giveConnExp_HCP(parc,op.tract,op.probe,weight,op.brainPart,op.nRem);

% 1. all genes: use exponential fit to correct for distance effect
CGEmatrix_uncorrected = corr(coexpData.parcelExpression(:,2:end)');
[CGEmatrixALL, FitCurve] = measure_correctDistance(CGEmatrix_uncorrected, coexpData.averageDistance, 'Correlated gene expression', 'exp', true);
CGEmatrixALL = tril(CGEmatrixALL);
CGEmatrixALL = CGEmatrixALL+CGEmatrixALL'; 

% 2. genes from enriched categories
% load selected GO categories:
GOcat = readtable('ermineJresultsRich_LRcortex_2.000000e-01perc_3.000000e-01_strength_k105_HCP_RNAseqprobes_NWDC_pearson_2020.csv');
GOcat = GOcat(GOcat.Pval_corr<0.05,:);
% import GEMMA annotation
GEMMA = importGEMMAfile('Generic_human_ncbiIds_noParents.an.txt');
genesSEL = cell(size(GOcat,1), 1);

for g=1:size(GOcat,1)
    GO = GOcat.GOcategory{g};
    GOf = [GO,'|'];
    GOe = ['|', GO];
    GOm = ['|', GO, '|'];
    % find all genes mentioned under this category
    G1 = find(contains(GEMMA.GOTerms, GOf));
    G2 = find(contains(GEMMA.GOTerms, GOe));
    G3 = find(contains(GEMMA.GOTerms, GOm));
    
    GIND = unique([G1; G2; G3]);
    
    genesSEL{g} = GEMMA.NCBIids(GIND);
    
end

genesALL = unique(vertcat(genesSEL{:}));

% select those genes from expression data
[~, INDexp] = intersect(coexpData.probeInformation.EntrezID, genesALL);

CGEmatrix_unENRICH = corr(coexpData.parcelExpression(:,INDexp+1)');
[CGEmatrixENRICH, FitCurve] = measure_correctDistance(CGEmatrix_unENRICH, coexpData.averageDistance, 'Correlated gene expression', 'exp', true);
CGEmatrixENRICH = tril(CGEmatrixENRICH);
CGEmatrixENRICH = CGEmatrixENRICH+CGEmatrixENRICH'; 

% Load BigBrain data
yregress = true;
[BBmean,microDist,distNorm, BBmpc, mp, BBskew] = getBBdata(parc, yregress);
BBmpc = BBmpc(1:size(CGEmatrixALL,1), 1:size(CGEmatrixALL,1)); 

save('data/modelling/CGEmatrices.mat', 'CGEmatrixALL', 'CGEmatrixENRICH', 'BBmpc')
end
