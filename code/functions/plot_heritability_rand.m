function plot_heritability_rand(n)

parcellation = 'HCP';
conWeight = 'FA';
op = selectCONmetrics(parcellation, conWeight);
whatFactors = {'Afactor'};
% A should be last - this matrix will be used for heritability plotting

plotOptions.colIn = [1 1 1];
plotOptions.colorOut = [82 82 82]/255; %[0,69,41]/255; %[.35 .35 .35]; %[4 90 141]/255; %[.45 .45 .45]; %[5 113 176]/255;
plotOptions.whatDistribution = 'histogram';


realTrajectory_save = cell(n,1); 

for k=1:length(whatFactors)
    figure('color','w');
    for jj=1:n
        
        [heritMatrix, nodeData, groupAdjlog, realTrajectory_save{jj}] = S3_compareHeritability_rand(jj, plotOptions);
        
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
        
        figureName = sprintf('makeFigures/rand/%s_curves_rand%d.png', whatFactors{k}, jj);
        print(gcf,figureName,'-dpng','-r600');
    end
    
% make a plot with empty bottom
plot_onlyHist(nodeData, plotOptions)

A = vertcat(realTrajectory_save{:}); 
myColors = GiveMeColors('RFPU'); 
sortK = sort(nodeData,'descend');
maxK = sortK(2); % Up to the second-highest k
kr = min(nodeData):maxK;
krAll = min(nodeData):max(nodeData);

for kk=1:3
    
    B = horzcat(A{:,kk});
    % plot R/F/P curves from all
    plot_distribution(kr,B', 'Color', myColors(kk,:));

end

set(gcf, 'Position', [500 500 750 500])

ylabel({'Mean h^2'})
xlabel('Node degree, k');
xlim([min(nodeData)-0.5,max(nodeData)+0.5]);
set(gca,'fontsize', 18);
ylim([0 0.05])

figureName = sprintf('makeFigures/rand/%s_curves_meanALL.png', whatFactors{k});
print(gcf,figureName,'-dpng','-r600');

end

