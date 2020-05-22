%--------------------------------------------------%
% Figure S2-S5: Heritability
%--------------------------------------------------%
function FigureS2_S5()

parcellation = 'HCP';
weight = 'FA';
op = selectCONmetrics(parcellation, weight);

plotOptions.colIn = [1 1 1];
plotOptions.colorOut = [82 82 82]/255; %[0,69,41]/255; %[.35 .35 .35]; %[4 90 141]/255; %[.45 .45 .45]; %[5 113 176]/255;
plotOptions.whatDistribution = 'histogram';

% S2 make a plot for ACTE model for all edges
heritMatrixACTE = S3_compareHeritability('HCP',op.tract,'Afactor', weight,op.densThreshold,op.cvMeasure, plotOptions, true);
ylim([0.2 0.6])
ylabel({'Mean h^2'})

figureName = sprintf('makeFigures/Afactorcurves_ACTEonly_%s_%d.png', 'HCP', 0.2);
print(gcf,figureName,'-dpng','-r600');


% S3 - heritability in distance ranges
% get heritability for best-fitting models
[heritMatrix, nodeData, groupAdjlog, mask] = S3_compareHeritability(parcellation,op.tract,'Afactor',weight,op.densThreshold,op.cvMeasure, plotOptions, false);
numThr = 4;
distMatr = giveConnDistance(parcellation, op.tract, groupAdjlog);

% bin data into distance bins
heritMatrixHalf = maskuHalf(heritMatrix);
heritMatrixHalf(groupAdjlog==0) = NaN;
distMatr = maskuHalf(distMatr);

% plot heritability for different distance bins as violin plots
[RvsF, FvsP, dataCell,~,f0] = plot_distanceViolin(heritMatrixHalf, distMatr, groupAdjlog, nodeData, op.khub, numThr, 'Heritability');
figureName = sprintf('makeFigures/heritability_distributions_distance_%s.png', parcellation);
print(f0,figureName,'-dpng','-r600');

% S4 - heritability for SC weight
whatFactors = {'Afactor', 'Efactor'};
weightSC = 'standard';

for k=1:length(whatFactors)
    
    [heritMatrixSC, nodeData, groupAdjlog, mask] = S3_compareHeritability('HCP',op.tract,whatFactors{k}, weightSC, op.densThreshold,op.cvMeasure, plotOptions, false);
    
    switch whatFactors{k}
        case 'Efactor'
            ylim([0.6 1])
            ylabel({'Mean e^2'})
            
        case 'Afactor'
            ylim([0 0.4])
            ylabel({'Mean h^2'})
    end
    
    figureName = sprintf('makeFigures/%s_curves_%s_%s_%d.png', whatFactors{k}, weightSC, parcellation, round(op.densThreshold*100));
    print(gcf,figureName,'-dpng','-r600');
    
end

% S5 - heritability for different parcellations and densities;
parcs = {'HCP','random500'};
tractography = 'iFOD2';
conWeight = 'FA';
densThresholds = [0.15, 0.2, 0.25; % for HCP parcellation
    0.05 0.1, 0.15]; % for random 500 parcellation
khubs = [90, 105, 120; % for HCP parcellation
    35, 70, 100]; % for random 500 parcellation
yVals = [[0.3, 0.7]; [0.3 0.7]];
cvMeasure = 'strength';

for pa=1:length(parcs)
    parc = parcs{pa};
    yVal = yVals(pa,:);
    
    for de=1:size(densThresholds,2)
        densThreshold = densThresholds(pa,de);
        khub = khubs(pa,de);
        
        [heritMatrix, nodeData, groupAdjlog, mask] = S3_compareHeritability(parc,tractography,'Afactor', conWeight,densThreshold,cvMeasure, plotOptions);
        
        figureName = sprintf('makeFigures/Afactor_curves_%s_%d.png', parc, densThreshold);
        ylim(yVal)
        ylabel({'Mean h^2'})
        
        print(gcf,figureName,'-dpng','-r600');
        
    end
end
end






