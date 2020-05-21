function [cellGenesUNIQUE,cellGroups] = getCellSpecificGenes(filename)

opts = detectImportOptions(filename);
opts = setvartype(opts,'string');

T = readtable(filename,opts);

% get a list for each group
cellGroups = unique(T.Var3(2:end)); 

for g=1:length(cellGroups)
    % for each group of genes, select unique gene names from the table and
    % remove missing entries
    indG = find(T.Var3==cellGroups(g)); 
    genesG = table2cell(T(indG,4:end)); 
    
     Y_temp = cellfun(@(M) M(:), genesG, 'Uniform', 0);
     Y = horzcat(Y_temp{:});
     R = unique(rmmissing(Y));
     
     cellGenes{g} = R; 
end

% from each list get genes that are unique to it and are not mentioned in
% others
for l=1:length(cellGroups)
    
    lind = zeros(length(cellGroups),1); 
    lind(l) = 1; 
    
    L_other = cellGenes(lind==0); 
    OTHER = horzcat(L_other{:});
    SELECT = cellGenes{lind==1}; 
    
    C = setdiff(SELECT,OTHER); 
    cellGenesUNIQUE{l} = C; 
    
end
end

