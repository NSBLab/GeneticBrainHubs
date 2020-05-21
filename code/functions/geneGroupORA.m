% ---------------------------------------------------------------------
% Find if genes that are outliers in non-cell specific list are
% ---------------------------------------------------------------------
function pORA = geneGroupORA(geneList, coexpData, listGenes, whatORA)

if nargin<4
    whatORA = 'NCELLgenes'; 
end

% 1. Find the outliers in group 10 (all the remaining genes)
switch whatORA
    case 'ALLgenes'
        ind10 = geneList.geneGroup~=0;
        y = quantile(geneList.geneScore(ind10),3);
        OUTall = find(geneList.geneScore>y(3));
        OUT10 = OUTall;
        lists = {'MET', 'HSE', 'HAR', 'scz', 'adhd', 'aut', 'bip', 'mdd', 'iq', 'Ex-Neuron', 'In-Neuron', 'OPC', 'astro', 'endo', 'micro', 'oligo'};
    case 'NCELLgenes'
        ind10 = geneList.geneGroup==10;
        y = quantile(geneList.geneScore(ind10),3);
        OUTall = find(geneList.geneScore>y(3)+1.5*iqr(geneList.geneScore(ind10)));
        %OUTall = find(geneList.geneScore>y(3));
        OUT10 = intersect(OUTall, find(geneList.geneGroup==10));
        lists = {'MET', 'HSE', 'HAR', 'scz', 'adhd', 'aut', 'bip', 'mdd', 'iq'};
end
% over-represented in pre-defined genes lists

% options specific for GWAS-based lists
brainParts = {'Cortex'}; %, 'Frontal_Cortex_BA9', 'Anterior_cingulate_cortex_BA24'};
numTissues = length(brainParts);

for l=1:length(lists)
    whatGeneSet = lists{l};
    if strcmp(whatGeneSet, 'MET')
        
        geneNames = importMouseGenes('mouse_metabolicGenes.txt');
        mgInd = cell(length(geneNames),1);
        for mg=1:length(geneNames)
            mgInd{mg} = find(strcmp(geneNames{mg},coexpData.probeInformation.GeneSymbol));
        end
        selectedGenes = cell2mat(mgInd);
    elseif strcmp(whatGeneSet, 'HSE')
        
        geneNames = {'BEND5', 'C1QL2', 'CACNA1E', 'COL24A1', 'COL6A1', 'CRYM', 'KCNC3', 'KCNH4', 'LGALS1', ...
            'MFGE8', 'NEFH', 'PRSS12', 'SCN3B', 'SCN4B', 'SNCG', 'SV2C', 'SYT2', 'TPBG', 'VAMP1'};
        mgInd = cell(length(geneNames),1);
        for mg=1:length(geneNames)
            mgInd{mg} = find(strcmp(geneNames{mg},coexpData.probeInformation.GeneSymbol));
        end
        selectedGenes = cell2mat(mgInd);
    elseif strcmp(whatGeneSet, 'HAR')
        % load gene list from
        HARgenes = importHARgenes('1-s2.0-S0092867416311692-mmc2.xlsx');
        kl=1;
        for g=2:length(HARgenes)
            if ischar(HARgenes{g})
                C = strsplit(HARgenes{g},',');
                for klp=1:length(C)
                    geneNames{kl} = C{klp};
                    kl=kl+1;
                end
            end
        end
        geneNames = unique(geneNames)';
        mgInd = cell(length(geneNames),1);
        for mg=1:length(geneNames)
            mgInd{mg} = find(strcmp(geneNames{mg},coexpData.probeInformation.GeneSymbol));
        end
        selectedGenes = cell2mat(mgInd);
        
        
    elseif strcmp(whatGeneSet, 'Ex-Neuron') || strcmp(whatGeneSet, 'In-Neuron') || strcmp(whatGeneSet, 'OPC')...
            || strcmp(whatGeneSet, 'astro') || strcmp(whatGeneSet, 'endo') || strcmp(whatGeneSet, 'micro') || strcmp(whatGeneSet, 'oligo')
        
        cellfind = @(string)(@(cell_contents)(strcmp(string,cell_contents)));
        groupID = find(cellfun(cellfind(whatGeneSet),listGenes(:,1)));
        [~, selectedGenes] = intersect(coexpData.probeInformation.EntrezID, geneList.geneID(geneList.geneGroup==groupID), 'stable');
        
        if strcmp(whatGeneSet, 'Ex-Neuron')
            whatGeneSet = 'ExNeuron';
        elseif strcmp(whatGeneSet, 'In-Neuron')
            whatGeneSet = 'InNeuron';
        end
        
    elseif strcmp(whatGeneSet, 'brainCombined')
        brainGenes = importBrainGenes('tissue_specificity_rna_cerebral_allCombined.txt');
        mgInd = cell(size(brainGenes,1),1);
        for mg=1:size(brainGenes,1)
            mgInd{mg} = find(strcmp(brainGenes{mg,1},probeInformation.GeneSymbol));
        end
        mgInd = mgInd(~cellfun('isempty',mgInd));
        selectGenes = cell2mat(mgInd);
        
    elseif strcmp(whatGeneSet, 'brainEnriched')
        brainGenes = importBrainGenes('tissue_specificity_rna_cerebral_groupEnriched.txt');
        mgInd = cell(size(brainGenes,1),1);
        for mg=1:size(brainGenes,1)
            mgInd{mg} = find(strcmp(brainGenes{mg,1},probeInformation.GeneSymbol));
        end
        mgInd = mgInd(~cellfun('isempty',mgInd));
        selectGenes = cell2mat(mgInd);
        
    elseif strcmp(whatGeneSet, 'brainEnhanced')
        brainGenes = importBrainGenes('tissue_specificity_rna_cerebral_tissueEnhanced.txt');
        mgInd = cell(size(brainGenes,1),1);
        for mg=1:size(brainGenes,1)
            mgInd{mg} = find(strcmp(brainGenes{mg,1},probeInformation.GeneSymbol));
        end
        mgInd = mgInd(~cellfun('isempty',mgInd));
        selectGenes = cell2mat(mgInd);
        
    elseif strcmp(whatGeneSet, 'brainCombinedHAR')
        
        % select brain genes
        brainGenes = importBrainGenes('tissue_specificity_rna_cerebral_allCombined.txt');
        mgInd = cell(size(brainGenes,1),1);
        for mg=1:size(brainGenes,1)
            mgInd{mg} = find(strcmp(brainGenes{mg,1},coexpData.probeInformation.GeneSymbol));
        end
        mgInd = mgInd(~cellfun('isempty',mgInd));
        selectGenesBRAIN = cell2mat(mgInd);
        
        % select HAR
        HARgenes = importHARgenes('1-s2.0-S0092867416311692-mmc2.xlsx');
        kl=1;
        for g=2:length(HARgenes)
            if ischar(HARgenes{g})
                C = strsplit(HARgenes{g},',');
                for klp=1:length(C)
                    geneNames{kl} = C{klp};
                    kl=kl+1;
                end
            end
        end
        
        geneNames = unique(geneNames)';
        mgInd = cell(length(geneNames),1);
        for mg=1:length(geneNames)
            mgInd{mg} = find(strcmp(geneNames{mg},coexpData.probeInformation.GeneSymbol));
        end
        selectedGenesHAR = cell2mat(mgInd);
        
        % find the overlap
        selectedGenes = intersect(selectGenesBRAIN, selectedGenesHAR);
        
    else
        results = selectGWASgenes(coexpData, whatGeneSet, brainParts, numTissues, 'oneList', 1);
        if sum(strcmp(fieldnames(results), 'geneID')) == 1
            [~, selectedGenes] = intersect(coexpData.probeInformation.EntrezID, results.geneID, 'stable');
        end
    end
    
    T = exist('selectedGenes');
    if T==1
        % X - number of genes that need to exceed
        % what is the actual number of genes from a list among outliers
        empOverlap = intersect(OUT10, selectedGenes);
        X = 0:length(empOverlap);
        % M - total population size
        M = length(find(ind10));
        % K - number of items with a desired characteristic (number of genes in a
        % specific list)
        K = length(selectedGenes);
        % N - number of genes selected
        N = length(OUT10);
        
        Y = hygepdf(X,M,K,N);
        pORA.(whatGeneSet) = 1-sum(Y);
    end
    clearvars selectedGenes
end
end
