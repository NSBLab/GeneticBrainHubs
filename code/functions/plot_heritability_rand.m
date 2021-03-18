function plot_heritability_rand(n)

parcellation = 'HCP';
conWeight = 'FA';
op = selectCONmetrics(parcellation, conWeight);
whatFactors = {'Afactor'};
% A should be last - this matrix will be used for heritability plotting

plotOptions.colIn = [1 1 1];
plotOptions.colorOut = [82 82 82]/255; %[0,69,41]/255; %[.35 .35 .35]; %[4 90 141]/255; %[.45 .45 .45]; %[5 113 176]/255;
plotOptions.whatDistribution = 'histogram';

for k=1:length(whatFactors)
    figure('color','w');
    for n=1:n
        
        [heritMatrix, nodeData, groupAdjlog] = S3_compareHeritability_rand(n, plotOptions);
        
        switch whatFactors{k}
            case 'Efactor'
                if strcmp(conWeight, 'FA')
                    ylim([0.3 0.7])
                else
                    ylim([0.5 1])
                end
                ylabel({'Mean e^2'})
                
            case 'Cfactor'
                ylim([0 0.1])
                ylabel({'Mean c^2'})
                
            case 'Afactor'
                if strcmp(conWeight, 'FA')
                    ylim([0 0.1])
                else
                    ylim([0 0.4])
                end
                ylabel({'Mean h^2'})
                
            case 'Tfactor'
                ylim([0 0.01])
                ylabel({'Mean phenotypic variance explained','by twin specific environmental factors'})
        end
        
        figureName = sprintf('makeFigures/rand/%s_curves_rand%d.png', whatFactors{k}, n);
        print(gcf,figureName,'-dpng','-r600');
    end
    
end

