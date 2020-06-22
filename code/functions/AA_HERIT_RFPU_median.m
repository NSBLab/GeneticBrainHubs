function [dataCell,pRF,statsRF,pFP,statsFP] = AA_HERIT_RFPU_median(Exp,AdjCon, nodeDeg, kHub)
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

% 
% colors = num2cell(rgb_colorMatrix, 2);
% extraParams.theColors = colors;
% extraParams.customSpot = '.';
% 
% nR = nansum(logical(dataCell{1}));
% nF = nansum(logical(dataCell{2}));
% nP = nansum(logical(dataCell{3}));
% 
% % f = figure('color', 'w');
% % f.Position = [1000,200,600,500];
% 
% JitteredParallelScatter(dataCell, true, true, false, extraParams)
% ylabel('Edge heritability', 'FontSize', 18);
% set(gca,'Xtick', [1 2 3], 'XTickLabel',...
% {sprintf('Rich (%d pairs)', nR), ...
% sprintf('Feeder (%d pairs)',nF),...
% sprintf('Peripheral (%d pairs)',nP)});
% set(gca,'fontsize',16)
% ylim([0 1])
% xtickangle(30)
% compare rich vs feeder links
[pRF,~,statsRF] = ranksum(richLinks,feederLinks, 'tail', 'right');
% compare feeder vs peripheral links
[pFP,~,statsFP] = ranksum(feederLinks,periphLinks, 'tail', 'right');

%t = stats.tstat;
end
