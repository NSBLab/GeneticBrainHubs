function [anovaOUT, geneList, fig] = GCCanova(coexpData, FitCurve, GrLC, nodeDegLC, whatGenes, kHub, listGenes)

indGENEsymb = vertcat(listGenes{:,4});
indGROUP = vertcat(listGenes{:,5});
nameGROUP = vertcat(listGenes{:,6});


% calculate gene contribution scores (GCC) and do ANOVA on them
switch whatGenes
    case 'ALL'
        % this calculates GCC using all available genes and compares
        % cell-specific with other genes in ANOVA
        indGENEsymbALL = 1:length(coexpData.probeInformation.DS);
        % calculate gene contribution score and compare them between link types
        GCC = giveGCC(coexpData.parcelExpression,indGENEsymbALL, FitCurve);
        geneScore = compareGCC(GCC, GrLC, 'Rich', nodeDegLC, kHub, true);
        geneID = coexpData.probeInformation.EntrezID(indGENEsymbALL);
        geneName = coexpData.probeInformation.GeneSymbol(indGENEsymbALL);
        
        % assign a number for genes in each category and a number that is not mentioned yet to other genes (e.g. 10)
        geneGroup = zeros(length(coexpData.probeInformation.DS),1);
        geneGroup(indGENEsymb) = indGROUP;
        geneGroup(geneGroup==0) = 10; 
        
        groupName = cell(length(coexpData.probeInformation.DS),1);
        groupName(indGENEsymb) = nameGROUP; 
        groupName(geneGroup==10) = {'all other genes'}; 
        
        geneList = table(geneID, geneName, geneScore,geneGroup, groupName);
        % parametric ANOVA (comparing means)
        [p,tbl,stats] = anova1(geneList.geneScore, geneList.geneGroup);
        groupNamesLabel = {'Ex-neuron'; 'In-neuron'; 'OPC'; 'astro'; 'endo'; 'micro'; 'oligo'; 'other'}; 
        xticklabels(groupNamesLabel);
        ylabel('Gene score (t-val)')
        
    case 'CELL'
        GCC = giveGCC(coexpData.parcelExpression,indGENEsymb, FitCurve);
        geneScore = compareGCC(GCC, GrLC, 'Rich', nodeDegLC, kHub, true);
        geneID = coexpData.probeInformation.EntrezID(indGENEsymb);
        
        geneList = table(geneID, geneScore,indGROUP,nameGROUP, indGENEsymb);
        % parametric ANOVA (comparing means)
        [p,tbl,stats] = anova1(geneList.geneScore, geneList.nameGROUP);
end

figure; 
[c,~,~,gnames] = multcompare(stats);

anovaOUT.p = p;
anovaOUT.tbl = tbl;
anovaOUT.stats = stats; 

% put values in a cell and make a violin plot
uG = [1:7,10];
allData = cell(length(uG),1);
groupNames = cell(length(uG),1);
k=1;
for i=uG
    ind = find(geneList.geneGroup==i);
    allData{k} = geneList.geneScore(ind);
    groupNames{k} = geneList.groupName{ind(1)};
    
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
ylabel('Gene contribution score (r>p)')

% change colour

groupsColor = [[33,113,181]/255; [8,69,148]/255; ...
   [161,217,155]/255; [116,196,118]/255; [65,171,93]/255; [35,139,69]/255; [0,109,44]/255; ...
   [254,153,41]/255]; 


violinLabels = cell(size(groupNames,1),1); 

for i=1:size(groupNames,1)
    violins(1,i).ViolinColor = groupsColor(i,:);
    violins(1,i).EdgeColor = groupsColor(i,:);
    violins(1,i).BoxColor = [.25 .25 .25];
    
    A = strcat(groupNamesLabel{i}, {' '}, {'('}, num2str(length(allData{i})),{')'}); 
    violinLabels{i} = A{1};

end

xticklabels(violinLabels); 
xtickangle(30)



end

