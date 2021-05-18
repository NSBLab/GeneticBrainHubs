function FigureS6_uncorrectedCGE()
%results presented in Figure S6 C

whatDWI = 'HCP';
weight = 'standard'; %for GenCog'standard'; 
parc = 'HCP';
op = selectCONmetrics(parc, weight); 

plotOptions.colIn = [1 1 1];
plotOptions.colorOut = [82 82 82]/255; %[0,69,41]/255; %[.35 .35 .35]; %[4 90 141]/255; %[.45 .45 .45]; %[5 113 176]/255; 
plotOptions.whatDistribution = 'histogram'; 

if strcmp(whatDWI, 'HCP')
    
    [coexpData, A, matrices, coordinates, avWeight] = giveConnExp_HCP(parc,op.tract,op.probe,weight,op.brainPart,op.nRem);
    
    
elseif strcmp(whatDWI, 'GenCog')

    [coexpData, A, matrices, coordinates, avWeight] = giveConnExp_GenCog(parc,op.tract,op.probe,weight,op.brainPart,op.nRem);
   
    
end

numLC = size(coexpData.averageCoexpression,1); 

[GrSC] = giveMeRichClub(matrices, coordinates, op.groupConn ,op.densThreshold, false, op.cvMeasure, op.consThr);
groupAdjlog = logical(GrSC); 
nodeData = degrees_und(groupAdjlog);


% select left hemisphere data for connectivity matrix, degree distribution
groupAdjlog = groupAdjlog(1:numLC, 1:numLC); 
nodeData = nodeData(1:numLC); 

% do not correct for distance effect
CGEmatrix = corr(coexpData.parcelExpression(:,2:end)');

CON = CGEmatrix(groupAdjlog==1); 
UNCON = CGEmatrix(groupAdjlog==0);
CON(isnan(CON)) = []; 
UNCON(isnan(UNCON)) = [];
[h,p,ci,stats] = ttest2(CON, UNCON, 'Vartype', 'unequal'); 

[dataCell,pRF,statsRF,pFP,statsFP, pRU,statsRU] = compare_CGE_RFPU(CGEmatrix,groupAdjlog, nodeData, op.khub, false); 
figureName = sprintf('makeFigures/CGEdistributions_uncorrected_%s_%d.png', parc,round(op.densThreshold*100));
print(gcf,figureName,'-dpng','-r300');

CGEmatrix(groupAdjlog==0) = NaN; 

% plot R/F/P lines for CGE
whatTail = 'right'; 
[~, data_export] = RichClubHuman(groupAdjlog,CGEmatrix, nodeData, whatTail, plotOptions.whatDistribution, plotOptions.colorOut, plotOptions.colIn); 
ylabel({'Mean correlated gene expression'})
set(gcf, 'Position', [500 500 750 550])
set(gca,'fontsize', 18);
ylim([0.05 0.45])

writetable(data_export,'data_export/source_data.xlsx','Sheet','Supplementary Figure6c','WriteVariableNames',true);
degree = nodeData';
region = (1:length(nodeData))';
node_degree = table(region,degree);
writetable(node_degree,'data_export/source_data.xlsx','Sheet','Supplementary Figure6c','Range','G:H','WriteVariableNames',true);


figureName = sprintf('makeFigures/CGEcurves_uncorrected_%s_%d.png', parc,round(op.densThreshold*100));
print(gcf,figureName,'-dpng','-r300');

end