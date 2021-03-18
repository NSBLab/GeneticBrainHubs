
load('100DS180scaledRobustSigmoidNSGRNAseqQC1Lcortex_ROI_NOdistCorrSurface.mat')

for i=1:5
    nodeData = parcelExpression(:,i+1); 
    plot_expressionSurface('HCP',nodeData,'outside', 'lh');
    figureName = sprintf('makeFigures/expression_gene%d.png', i);
    print(gcf,figureName,'-dpng','-r300');
end