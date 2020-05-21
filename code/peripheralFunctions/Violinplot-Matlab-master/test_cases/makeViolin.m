% make violin plots for modules
% put data to structure; 
close all; clear all; 
plotWhat = 'heritability'; 
load('matrixHERITABILITY.mat')

switch plotWhat
    case 'CGE'
        load('matrixCGE.mat')
        load('modulesCGE.mat')
        nodeData = nodeData(1:180); 
        adj = groupAdjlog(1:180, 1:180); 
        matrix = CGEmatrix; 
        yName = 'Correlated gene expression'; 
        yRange = [-1 1]; 
    case 'heritability'
        adj = groupAdjlog; 
        load('modulesHERITABILITY.mat')
        yName = 'Heritability'; 
        yRange = [0 1]; 
        heritMatrix(heritMatrix==0) = NaN; 
        matrix = heritMatrix; 
end

modulesColor = vertcat( [254,224,144]/255,  [253,174,97]/255,[77,77,77]/255, ...
                     [178,24,43]/255,[135,135,135]/255, [244,165,130]/255, ...
                     [146,197,222]/255,[67,147,195]/255, [33,102,172]/255, ...
                     [27,120,55]/255, [214,96,77]/255, [224,224,224]/255); 
                 
for i=1:12
    switch plotWhat
        case 'heritability'
            moduleVals = allData{i};
            moduleVals(moduleVals==0) = NaN;
            S.(netNamesAll{i}) = moduleVals;
        case 'CGE'
            S.(netNamesAll{i}) = allData{i};
    end
    
end

set(gca,'FontSize',18)
violins = violinplot(S);  
set(gcf,'color','w');
ylabel(yName)
ylim(yRange)
% change colors
for i=1:size(modulesColor,1)
violins(1,i).ViolinColor = modulesColor(i,:); 
end


%% make violin plots for inter-moduler and intra-modular rich links
[MeasureINT] = compareRichModules(matrix, adj, nodeData, 105, plotWhat, yRange); 
close all; 

netNamesAll = {'INTRArich'; 'INTRAnonrich'; 'INTERrich'; 'INTERnonrich'};
c = vertcat([227 26 28]/255, [158 202 225]/255, [227 26 28]/255, [158 202 225]/255); 

for j=1:length(MeasureINT)
    J.(netNamesAll{j}) = MeasureINT{j}; 
end

set(gca,'FontSize',18)
violins = violinplot(J);
set(gcf,'color','w');
ylabel(yName)
% change colors
for i=1:size(c,1)
violins(1,i).ViolinColor = c(i,:); 
end
xticklabels({'Rich', 'Non-rich', 'Rich', 'Non-rich'})
ylim(yRange)

%% make plot for GCC groups
load('cellSpecificGCC.mat')
cellGroups = {'ExNeuron', 'InNeuron', 'OPC', 'astro','endo','micro','oligo','other'}; 
for k=1:length(groupNames)
    GCC.(cellGroups{k}) = allData{k};
end
set(gca,'FontSize',18)
violins = violinplot(GCC);
set(gcf,'color','w');
xticklabels(groupNames)
ylabel('Gene contribution score (rich>peripheral)')

%% plot modules on surface
load('cortex_parcel_network_assignments.mat')

[f] = plot_measureONsurface('HCP', netassignments, 'outside', modulesColor)
[f] = plot_measureONsurface('HCP', netassignments, 'inside', modulesColor)
