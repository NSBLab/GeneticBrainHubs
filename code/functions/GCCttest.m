function [T, geneList, fig] = GCCttest(coexpData, FitCurve, GrLC, nodeDegLC, kHub, listGenes)

groupName = vertcat(listGenes{:,1});
indGENEsymb = vertcat(listGenes{:,4});
indGROUP = vertcat(listGenes{:,5});


indGENEsymbALL = 1:length(coexpData.probeInformation.DS);
% calculate gene contribution score and compare them between link types
GCC = giveGCC(coexpData.parcelExpression,indGENEsymbALL, FitCurve);
geneScore = compareGCC(GCC, GrLC, 'Rich', nodeDegLC, kHub, true);
        
% for each gene group, compare cell-specific vs all other genes
T = table; 
geneGroup = cell(length(geneScore),1); 
groupID = zeros(length(geneScore),1); 

for g=1:length(unique(indGROUP))
    
    iG = indGENEsymb(indGROUP==g); 
    cellGROUP = geneScore(iG); 
    iO = setdiff(indGENEsymbALL, iG); 
    otherGROUP = geneScore(iO); 
    [~, pG,~,stats] = ttest2(cellGROUP,otherGROUP, 'Vartype','unequal', 'Tail', 'right'); 
    nG = groupName(g);
    
    geneGroup(iG) = {nG}; 
    groupID(iG) = g; 
    T.group(g,1) = nG; 
    T.p(g,1) = pG; 
    T.tval(g,1) = stats.tstat; 

end


otherIND = find(cellfun('isempty', geneGroup)); 
geneGroup(otherIND) = {'other'}; 

otherIND2 = find(groupID==0); 
groupID(otherIND2) = size(listGenes,1)+1; 

geneName = coexpData.probeInformation.GeneSymbol; 
geneList = table(geneName, geneScore, geneGroup, groupID);

% put values in a cell and make a violin plot
uG = 1:length(unique(groupID)); 
allData = cell(length(uG),1);
groupNames = cell(length(uG),1);
k=1;
for i=uG
    ind = find(geneList.groupID==i);
    allData{k} = geneList.geneScore(ind);
    groupNames{k} = geneList.geneGroup{ind(1)};
    
    isD = contains(groupNames{k},'-');
    isS = contains(groupNames{k},' ');
    
    if isD
        groupNames{k} = erase(groupNames{k},'-');
    end
    
    if isS
        groupNames{k} = erase(groupNames{k},' ');
    end
        
    
    N = char(groupNames{k});
    S.(N) = allData{k};
    k=k+1;
    
end

% make a violin plot
fig = figure('color','white');
set(gcf, 'Position', [10 10 1200 500]);

set(gca,'FontSize',18)
violins = violinplot(S);
ylabel('GCS_{t-stat}')

% change colour

groupsColor = [[33,113,181]/255; [8,69,148]/255; ...
   [161,217,155]/255; [116,196,118]/255; [65,171,93]/255; [35,139,69]/255; [0,109,44]/255; ...
   [252,146,114]/255; [239,59,44]/255; [252,187,161]/255; [203,24,29]/255; [251,106,74]/255; [254,224,210]/255]; 
   

violinLabels = cell(size(groupNames,1),1); 

for i=1:size(groupNames,1)
    violins(1,i).ViolinColor = groupsColor(i,:);
    violins(1,i).EdgeColor = groupsColor(i,:);
    violins(1,i).BoxColor = [.25 .25 .25];
    
    A = strcat(groupNames{i}, {' '}, {'('}, num2str(length(allData{i})),{')'}); 
    violinLabels{i} = A{1};

end

xticklabels(violinLabels); 
xtickangle(30)

end

