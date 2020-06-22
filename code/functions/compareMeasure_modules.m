function [res,netNamesAll,allData, measureModRICH] = compareMeasure_modules(pairwiseMeasure,netassignments, groupAdjlog, nodeDeg, hubThr)
%---------------------------------------------------------------------------
% pairwiseMeasure, groupAdjlog, nodeDeg, hubThr, measureName, yVals

% This script compares heritability between modules
% the modular assignment defined based on functional networks
% https://github.com/ColeLab/ColeAnticevicNetPartition
%
% /Users/Aurina/GoogleDrive/Genetics_connectome/Heritability/data/connectomes/modules/cortex_parcel_network_assignments.mat
% defines modules.
%---------------------------------------------------------------------------
% pairwiseMeasure, groupAdjlog, nodeDeg, hubThr, measureName, yVals

% This script compares heritability between modules
% the modular assignment defined based on functional networks
% https://github.com/ColeLab/ColeAnticevicNetPartition
%
% /Users/Aurina/GoogleDrive/Genetics_connectome/Heritability/data/connectomes/modules/cortex_parcel_network_assignments.mat
% defines modules.

% ------
% INPUTS
% ------
% pairwiseMeasure   - pairwise measure in a region x region format
% groupAdjlog       - pairwise connectivity matrix in a region x region format
% nodeDeg           - node degree for each region in the connectivity matrix
% hubThr            - selected threshold for defininf hubs
% measureName       - pairwise measure name
% yVals             - range of y values to be used in the plot

%
% -------
% OUTPUTS
% -------
% res               - structure containing results: p and z value from a ranksum test comparing RICH links to other link groups. 
% netNamesAll       - cell containing module names
% allData           - cell containing pairwise measure for each module

%---------------------------------------------------------------------------

%---------------------------------------------------------------------------
% Load modular assignment
%---------------------------------------------------------------------------

netNamesAll = {'VIS1';'VIS2';'SM';'CO';'DAN';'LAN';'FPN';'AUD';'DMN';'PM';'VM';'OA'; 'RICH'; 'FEEDER'; 'PERIPH'; 'INTRA'; 'INTER'};
if size(pairwiseMeasure,1)==180
    netassignments = netassignments(1:180);
end
% keep only data for connected links and exclude zero heritabiltiy values
pairwiseMeasure(groupAdjlog==0) = NaN;
%---------------------------------------------------------------------------
% Extract heritability values for each link type: rich/feeder/peripheral
%---------------------------------------------------------------------------
[rgb_colorMatrix] = GiveMeColors('RFPU');
[measureRFP] = AA_HERIT_RFPU(pairwiseMeasure,groupAdjlog, nodeDeg, hubThr);

%---------------------------------------------------------------------------
% Extract heritability values for each module
%---------------------------------------------------------------------------
numModules = length(unique(netassignments));
measureMod = cell(numModules,1);

for m=1:numModules
    modVals = pairwiseMeasure(netassignments==m, netassignments==m);
    measureMod{m} = modVals(~isnan(modVals(:)));
end

% make a separate cell to contain rich links within each module so they can be plotted separately. 
measureModRICH = cell(numModules,1);
isHub = nodeDeg>hubThr; 
isRich = isHub.*isHub'.*groupAdjlog; 
isRich(isRich==0) = NaN; 
richLinks = pairwiseMeasure.*isRich;

for m=1:numModules
    modVals = richLinks(netassignments==m, netassignments==m);
    measureModRICH{m} = modVals(~isnan(modVals(:)));
end


%---------------------------------------------------------------------------
% Extract heritability values for inter-modular and intra-modular links
%---------------------------------------------------------------------------

intraModule = repmat(netassignments,[1,length(netassignments)])==repmat(netassignments',[length(netassignments),1]);
interModule = ~intraModule;

intraMeasure = pairwiseMeasure(intraModule==1);
interMeasure = pairwiseMeasure(interModule==1);

MeasureINT{1} = intraMeasure(~isnan(intraMeasure));
MeasureINT{2} = interMeasure(~isnan(interMeasure));

allData = vertcat(measureMod,measureRFP,MeasureINT');

% JitteredParallelScatter(allData, true, true, true)
% ylabel(sprintf('Edge %s',measureName), 'FontSize', 18);
% xlabel('Modules', 'FontSize', 18);
% set(gca,'Xtick', [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17], 'XTickLabel',netNamesAll);
% set(gca,'fontsize',16)
% ylim(yVals)

%---------------------------------------------------------------------------
% Compare heritability of rich links with each individual category using
% non-parametric test
%---------------------------------------------------------------------------
measureCAT = vertcat(measureMod,measureRFP,MeasureINT');
Rind = find(strcmp(netNamesAll, 'RICH'));

for c=1:length(measureCAT)
    if c~=Rind
        if ~isnan(measureCAT{Rind})
            [p,~,stats] = ranksum(measureCAT{Rind},measureCAT{c}, 'tail', 'right');
            res.(netNamesAll{c}).p = p;
            res.(netNamesAll{c}).z = stats.zval;
            fprintf('%s vs %s: p=%d, z=%d\n', netNamesAll{Rind}, netNamesAll{c}, round(p,4), round(stats.zval,4))
        end
    end
end

end
