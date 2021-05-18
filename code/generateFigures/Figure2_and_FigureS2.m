%--------------------------------------------------%
% Figure 2
%--------------------------------------------------%
function Figure2_and_FigureS2()

parcellation = 'HCP';
conWeight = 'FA';
op = selectCONmetrics(parcellation, conWeight);
whatFactors = {'Efactor', 'Afactor'};
% A should be last - this matrix will be used for heritability plotting

plotOptions.colIn = [1 1 1];
plotOptions.colorOut = [82 82 82]/255; %[0,69,41]/255; %[.35 .35 .35]; %[4 90 141]/255; %[.45 .45 .45]; %[5 113 176]/255;
plotOptions.whatDistribution = 'histogram';

for k=1:length(whatFactors)
    
    [heritMatrix, nodeData, groupAdjlog, mask, data_export] = S3_compareHeritability(parcellation,op.tract,whatFactors{k},op.weight,op.densThreshold,op.cvMeasure, plotOptions, false, 1000);
    
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
    
    % save data export to excel
    switch whatFactors{k}
        case 'Efactor'
    writetable(data_export,'data_export/source_data.xlsx','Sheet','Figure2d','WriteVariableNames',true);
        case 'Afactor'
    writetable(data_export,'data_export/source_data.xlsx','Sheet','Figure2c','WriteVariableNames',true);
    end
    
end
degree = nodeData';
region = (1:length(nodeData))';  
node_degree = table(region,degree); 
writetable(node_degree,'data_export/source_data.xlsx','Sheet','Figure2b','WriteVariableNames',true);

heritMatrixHalf = maskuHalf(heritMatrix);
heritMatrixHalf(groupAdjlog==0) = NaN;

% instead of ploting all connections - plot TOP 1000 heritability values
% and LOWEST 1000 heritability values; 
% highest heritability
f=figure('color','white');
set(gcf, 'Position', [10 10 700 1000]);
maskFULL = mask+mask';
HIher = double(maskFULL==2); 
HImap = [[215,48,31]/255; [215,48,31]/255]; 
plot_edges_brain(HIher,nodeData.^2,HImap); 
figureName = 'makeFigures/heritability_highest_brain_top.png';
print(f,figureName,'-dpng','-r300');

view([90 0]); camlight('right')
figureName = 'makeFigures/heritability_highest_brain_side.png';
print(f,figureName,'-dpng','-r300');

f=figure('color','white');
set(gcf, 'Position', [10 10 700 1000]);
LOher = double(maskFULL==1); 
LOmap = [[33,113,181]/255; [33,113,181]/255]; 
plot_edges_brain(LOher,nodeData.^2,LOmap); 
figureName = 'makeFigures/heritability_lowest_brain_top.png';
print(f,figureName,'-dpng','-r300');

view([90 0]); camlight('right')
figureName = 'makeFigures/heritability_lowest_brain_side.png';
print(f,figureName,'-dpng','-r300');


% plot violin distributions for modules add brain representation
load('data/modules/cortex_parcel_network_assignments.mat')
[resMod,netNamesAll,allData, measureModRICH] = compareMeasure_modules(heritMatrixHalf, netassignments, groupAdjlog, nodeData, op.khub);
[resINTERINTRA,f,S] = plot_moduleViolin(allData, measureModRICH, netNamesAll, netassignments, heritMatrixHalf, nodeData, groupAdjlog, 'heritability', op.khub);

% prepare data for export
num_points = cellfun(@length,struct2cell(S)); 
num_violins = length(fieldnames(S)); 
names_fields = fieldnames(S); 

S_exp = nan(max(num_points),num_violins);
for pp=1:num_violins
    S_exp(1:num_points(pp),pp) = S.(names_fields{pp});
end

S_export = array2table(S_exp,'VariableNames',names_fields); 
writetable(S_export,'data_export/source_data.xlsx','Sheet','Figure2g-i','WriteVariableNames',true);

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
[F, Hp, k_range] = plot_proportionHubsModule(netassignments, nodeData,plotOptions);
figureName = sprintf('makeFigures/proportion_HubsModules.png');
print(F,figureName,'-dpng','-r600');

% prepare data for exort
degree_threshold = [k_range; Hp];
net_names = [{''};names_fields(1:12)]; 
Hp_export = table(net_names,degree_threshold); 
writetable(Hp_export,'data_export/source_data.xlsx','Sheet','Figure2f','WriteVariableNames',true);

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
        
        plot_hubsSurface_SO(parcellation,ds,ts, side, hemi);
        
        figureName = sprintf('makeFigures/hubsSurface_%s_%d_%s_%s.png', parcellation, round(op.densThreshold*100), side, hemi);
        print(gcf,figureName,'-dpng','-r300');
    end
end
end



