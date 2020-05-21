function [dataCell,pRF,statsRF,pFP,statsFP] = AA_HERIT_RFPU(Exp,AdjCon, nodeDeg, kHub)
%
% connected vs unconnected
% Rexp - coexpression matrix
% Adj - adjecency matrix

%% output
% y - mean values for connected and unconnected links
% e - SD for connected and unconnected links
% dataCell{1,1} - distribution of connected values
% dataCell{2,1} - - distribution of unconnected values
% t - ttest2 t value
% p - ttest2 p value


mask = zeros(size(AdjCon));
isHub = nodeDeg>kHub;

mask(isHub, isHub) = 1; % rich
mask(isHub, ~isHub) = 2; % feeder
mask(~isHub, isHub) = 2; % feeder
mask(~isHub, ~isHub) = 3; % peripheral

AdjCon = double(AdjCon); 
AdjCon(AdjCon==0) = NaN;
conLinks = maskuHalf(AdjCon.*Exp);


richLinks = conLinks(mask==1); 
richLinks(isnan(richLinks)) = []; 

feederLinks = conLinks(mask==2); 
feederLinks(isnan(feederLinks)) = [];

periphLinks = conLinks(mask==3); 
periphLinks(isnan(periphLinks)) = [];


dataCell{1,1} = richLinks;
dataCell{2,1} = feederLinks;
dataCell{3,1} = periphLinks;

% compare rich vs feeder links
[hRF,pRF,ciRF,statsRF] = ttest2(richLinks,feederLinks, 'Vartype','unequal', 'tail', 'right');
% compare feeder vs peripheral links
[hFP,pFP,ciFP,statsFP] = ttest2(feederLinks,periphLinks, 'Vartype','unequal', 'tail', 'right');

%t = stats.tstat;
end
