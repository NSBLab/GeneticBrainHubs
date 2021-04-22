function plot_genetic_variance()

parcellation = 'HCP';
conWeight = 'FA';
op = selectCONmetrics(parcellation, conWeight);
whatFactors = {'Avariance'};
% A should be last - this matrix will be used for heritability plotting

plotOptions.colIn = [1 1 1];
plotOptions.colorOut = [82 82 82]/255; %[0,69,41]/255; %[.35 .35 .35]; %[4 90 141]/255; %[.45 .45 .45]; %[5 113 176]/255;
plotOptions.whatDistribution = 'histogram';

[heritMatrix, nodeData, groupAdjlog, mask] = S3_compareHeritability_var(parcellation,op.tract,whatFactors{1},op.weight,op.densThreshold,op.cvMeasure, plotOptions, false);

print(gcf,'makeFigures/genetic_variance.png','-dpng','-r600');

end