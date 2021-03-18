%--------------------------------------------------%
% Figure 2
%--------------------------------------------------%
function plot_tstat_ratio()

parcellation = 'HCP';
conWeight = 'FA';
whatFactors = {'Afactor'};
numPerm = 5000; 
% A should be last - this matrix will be used for heritability plotting

plotOptions.colIn = [1 1 1];
plotOptions.colorOut = [82 82 82]/255; %[0,69,41]/255; %[.35 .35 .35]; %[4 90 141]/255; %[.45 .45 .45]; %[5 113 176]/255;
plotOptions.whatDistribution = 'histogram';

for k=1:length(whatFactors)
    
    [heritMatrix, nodeData, groupAdjlog] = S3_compareHeritability_perm(plotOptions, whatFactors{k}, numPerm); 
    
    figureName = sprintf('makeFigures/log_tstat_ratio_nperm%d.png', numPerm);
    print(gcf,figureName,'-dpng','-r600');
    
    
end
end



