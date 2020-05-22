%--------------------------------------------------%
% Figure 2
%--------------------------------------------------%
function Figure2()

parcellation = 'HCP';
conWeight = 'FA';
op = selectCONmetrics(parcellation, conWeight);
whatFactors = {'Efactor', 'Afactor'};
% A should be last - this matrix will be used for heritability plotting

plotOptions.colIn = [1 1 1];
plotOptions.colorOut = [82 82 82]/255; %[0,69,41]/255; %[.35 .35 .35]; %[4 90 141]/255; %[.45 .45 .45]; %[5 113 176]/255;
plotOptions.whatDistribution = 'histogram';

for k=1:length(whatFactors)
    
    [heritMatrix, nodeData, groupAdjlog, mask] = S3_compareHeritability(parcellation,op.tract,whatFactors{k},op.weight,op.densThreshold,op.cvMeasure, plotOptions, false);
    
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
                ylim([0.3 0.7])
            else
                ylim([0 0.4])
            end
            ylabel({'Mean h^2'})
            
        case 'Tfactor'
            ylim([0 0.01])
            ylabel({'Mean phenotypic variance explained','by twin specific environmental factors'})
    end
    
    figureName = sprintf('makeFigures/%s_curves_%s_%s_%d.png', whatFactors{k}, conWeight, parcellation, round(op.densThreshold*100));
    print(gcf,figureName,'-dpng','-r600');
end

heritMatrixHalf = maskuHalf(heritMatrix);
heritMatrixHalf(groupAdjlog==0) = NaN;

% plot violin distributions for modules add brain representation
load('data/modules/cortex_parcel_network_assignments.mat')
[resMod,netNamesAll,allData, measureModRICH] = compareMeasure_modules(heritMatrixHalf, netassignments, groupAdjlog, nodeData, op.khub);
[resINTERINTRA,f] = plot_moduleViolin(allData, measureModRICH, netNamesAll, netassignments, heritMatrixHalf, nodeData, groupAdjlog, 'heritability', op.khub);

figureName = 'makeFigures/heritability_modules_distributions.png';
print(f,figureName,'-dpng','-r600');

[f1, f2, f3, f4] = plot_modulesSurface(netassignments);

figureName = sprintf('makeFigures/MOD_brainL_outside_%s.png', parcellation);
print(f1,figureName,'-dpng','-r300');

figureName = sprintf('makeFigures/MOD_brainR_inside_%s.png', parcellation);
print(f2,figureName,'-dpng','-r300');

figureName = sprintf('makeFigures/MOD_brainL_inside_%s.png', parcellation);
print(f3,figureName,'-dpng','-r300');

figureName = sprintf('makeFigures/MOD_brainR_outside_%s.png', parcellation);
print(f4,figureName,'-dpng','-r300');

% plot the proportion of hubs in each module as a function of degree
F = plot_proportionHubsModule(netassignments, nodeData,plotOptions);
figureName = sprintf('makeFigures/proportion_HubsModules.png');
print(F,figureName,'-dpng','-r600');


% plot degree on the cortical surface for 3 hub thresholds
ts = [145,125,105]; % hub thresholds, degree at which regions are labeled hubs;
sides = {'inside'; 'outside'};
hemis = {'rh'; 'lh'};
for s=1:2
    side = sides{s};
    
    for h=1:2
        hemi = hemis{h};
        
        if strcmp(hemi, 'lh')
            ds = nodeData(1:length(nodeData)/2);
        elseif strcmp(hemi, 'rh')
            ds = nodeData(length(nodeData)/2+1:length(nodeData));
        end
        
        plot_hubGroupsSurface(parcellation,ds,ts, side, hemi);
        
        figureName = sprintf('makeFigures/hubsSurface_%s_%d_%s_%s.png', parcellation, round(op.densThreshold*100), side, hemi);
        print(gcf,figureName,'-dpng','-r300');
    end
end
end



