%--------------------------------------------------%
% Figure 2
%--------------------------------------------------%
function plot_herit_permPval()

whatFactors = {'Afactor'};
numPerm = 5000; 
% A should be last - this matrix will be used for heritability plotting

plotOptions.colIn = [1 1 1];
plotOptions.colorOut = [82 82 82]/255; %[0,69,41]/255; %[.35 .35 .35]; %[4 90 141]/255; %[.45 .45 .45]; %[5 113 176]/255;
plotOptions.whatDistribution = 'histogram';

for k=1:length(whatFactors)
    
    S3_compareHeritability_permPval(plotOptions, whatFactors{k}, numPerm); 
    print(gcf,'makeFigures/heritability_curves_permPval.png','-dpng','-r600');
    
    
end
end



