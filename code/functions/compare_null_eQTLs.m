% compare real data overlap with GWAS genes copared to null
function R = compare_null_eQTLs(whatDATA)

if nargin < 1
    whatDATA = 'HCP';
end

% load results for real data
load('POSvsNEG_comparison.mat')
switch whatDATA
    case 'HCP'
        pORA_real = HCP_pORA_eQTL_POS;
    case 'GenCog'
        pORA_real = GenCog_pORA_eQTL_POS;
end
D = fields(pORA_real);

% get all file names in the permutation file
fileList = dir(sprintf('data/reeqtls/permutations/%s*', whatDATA));

Dall = nan(length(fileList), length(D));

% for each GWAS list make a distribution
for i = 1:length(fileList)
    
    load(fileList(i).name)
    switch whatDATA
        case 'HCP'
            pORA_null = HCP_pORA_eQTL;
        case 'GenCog'
            pORA_null = GenCog_pORA_eQTL;
    end
    
    % select nulls only for relevant cases
    [Dt,Di] = intersect(D, fields(pORA_null), 'stable');
    
    for j = 1:length(Dt)
        
        Dall(i,Di(j)) = pORA_null.(Dt{j}).Noverlap;
    end
    
end

% calculate p-value for each individual disorder
pVal = nan(length(D),1);
pVal_FDR = nan(length(D),1);

for k = 1:length(D)
    nR = pORA_real.(D{k}).Noverlap;
    
    % select only nulls with non-nan values
    Dsel = Dall(:,k);
    Dsel(isnan( Dsel)) = [];

    if ~isempty(Dsel)
        pVal(k) = mean(Dsel>=nR);
    end
    clearvars 'Dsel'
    
    % now FDR-corrected - select max from each null
    DFDR = max(Dall,[], 2);
    % select only non-nans just in case there are some
    DFDR(isnan(DFDR)) = [];
    
    if ~isempty(DFDR)
        pVal_FDR(k) = mean(DFDR>=nR);
    end
    clearvars 'DFDR'
end

R = table(D, pVal, pVal_FDR);

end
