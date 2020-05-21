% for each cell-specific list of genes find ones in the expression data
function listGenes = getCellExpression(cellGenesUNIQUE,cellGroups,coexpData)
k=1;
for cl=1:length(cellGroups)
    
    [~,~,indEXP] = intersect(cellGenesUNIQUE{cl}, coexpData.probeInformation.GeneSymbol, 'stable');
    fprintf('Group %s: %d/%d genes found\n', cellGroups(cl), length(indEXP), length(cellGenesUNIQUE{cl}))
    % combine expression values of those genes into a single matrix and
    % make a new list of genes in each category
    if ~isempty(indEXP)
        
        groupIND = zeros(length(indEXP),1);
        groupIND(groupIND==0) = k;
        
        groupName = cell(length(indEXP),1);
        groupName(cellfun('isempty',groupName)) = {cellGroups(cl)};
        
        listGenes{k,1} = cellGroups(cl);
        listGenes{k,2} = coexpData.probeInformation.GeneSymbol(indEXP);
        listGenes{k,3} = coexpData.parcelExpression(:,indEXP+1);
        listGenes{k,4} = indEXP;
        listGenes{k,5} = groupIND;
        listGenes{k,6} = groupName;
        k=k+1;
    end
    
end
end