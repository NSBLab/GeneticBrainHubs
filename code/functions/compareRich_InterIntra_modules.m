function [MeasureINT, res, MeasureINTnames] = compareRich_InterIntra_modules(pairwiseMeasure, netassignments, groupAdjlog, nodeDeg, hubThr)
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
if size(pairwiseMeasure,1)==180
    netassignments = netassignments(1:180);
end
% keep only data for connected links
pairwiseMeasure(groupAdjlog==0) = NaN;

%---------------------------------------------------------------------------
% Extract heritability values for inter-modular and intra-modular links
%---------------------------------------------------------------------------

intraModule = repmat(netassignments,[1,length(netassignments)])==repmat(netassignments',[length(netassignments),1]);
interModule = ~intraModule;
% get intra-modular and rich
isHub = nodeDeg>hubThr; 
hubMask = isHub.*isHub'; 
periphbMask = ~isHub.*~isHub; 
RP = hubMask | periphbMask; 

feedMask = ~RP; 
%~hubMask; 


intraRICHmask = intraModule.*hubMask; 
intraPERIPHmask = intraModule.*periphbMask; 
intraFEEDmask = intraModule.*feedMask; 

interRICHmask = interModule.*hubMask; 
interPERIPHmask = interModule.*periphbMask; 
interFEEDmask = interModule.*feedMask; 

% get intra-modular and non-rich

intraRICH = pairwiseMeasure(intraRICHmask==1);
interRICH = pairwiseMeasure(interRICHmask==1);

intraPERIPH = pairwiseMeasure(intraPERIPHmask==1);
interPERIPH = pairwiseMeasure(interPERIPHmask==1);

intraFEED = pairwiseMeasure(intraFEEDmask==1);
interFEED = pairwiseMeasure(interFEEDmask==1);


MeasureINT{1} = intraRICH(~isnan(intraRICH));
MeasureINT{2} = intraFEED(~isnan(intraFEED));
MeasureINT{3} = intraPERIPH(~isnan(intraPERIPH));

MeasureINT{4} = interRICH(~isnan(interRICH));
MeasureINT{5} = interFEED(~isnan(interFEED));
MeasureINT{6} = interPERIPH(~isnan(interPERIPH));

MeasureINTnames = {'intraRich', 'intraFeeder', 'intraPeripheral', 'interRich', 'interFeeder', 'interPeripheral'}; 

% extraParams.theColors{1,:} = [0.8 0.2 0.2]; 
% extraParams.theColors{2,:} = [0.11 0.22 0.73]; 
% extraParams.theColors{3,:} = [0.8 0.2 0.2]; 
% extraParams.theColors{4,:} = [0.11 0.22 0.73]; 

% use t-test here
[~,p,~,stats] = ttest2(MeasureINT{1},MeasureINT{2}, 'Vartype','unequal','tail', 'right');
res.INTRA.pRF = p;
res.INTRA.tRF = stats.tstat;

[~,p,~,stats] = ttest2(MeasureINT{1},MeasureINT{3}, 'Vartype','unequal','tail', 'right');
res.INTRA.pRP = p;
res.INTRA.tRP = stats.tstat;

[~,p,~,stats] = ttest2(MeasureINT{4},MeasureINT{5}, 'Vartype','unequal','tail', 'right');
res.INTER.pRF = p;
res.INTER.tRF = stats.tstat;

[~,p,~,stats] = ttest2(MeasureINT{4},MeasureINT{6}, 'Vartype','unequal','tail', 'right');
res.INTER.pRP = p;
res.INTER.tRP = stats.tstat;

end
