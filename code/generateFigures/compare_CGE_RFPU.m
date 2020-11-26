function [dataCell,pRF,statsRF,pFP,statsFP, pRU,statsRU] = compare_CGE_RFPU(Exp,AdjCon, nodeDeg, kHub, isCorrected)
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
AdjUncon = double(~AdjCon); 
AdjCon(AdjCon==0) = NaN;
AdjUncon(AdjUncon==0) = NaN;

conLinks = maskuHalf(AdjCon.*Exp);
unconLinks = maskuHalf(AdjUncon.*Exp);
unconLinks(isnan(unconLinks)) = []; 

richLinks = conLinks(mask==1); 
richLinks(isnan(richLinks)) = []; 

feederLinks = conLinks(mask==2); 
feederLinks(isnan(feederLinks)) = [];

periphLinks = conLinks(mask==3); 
periphLinks(isnan(periphLinks)) = [];


dataCell{1,1} = richLinks;
dataCell{2,1} = feederLinks;
dataCell{3,1} = periphLinks;
dataCell{4,1} = unconLinks; 

% compare rich vs feeder links
[~,pRF,~,statsRF] = ttest2(richLinks,feederLinks, 'Vartype','unequal', 'tail', 'right');
% compare feeder vs peripheral links
[~,pFP,~,statsFP] = ttest2(feederLinks,periphLinks, 'Vartype','unequal', 'tail', 'right');

% compare rich vs unconnected
[~,pRU,~,statsRU] = ttest2(richLinks,unconLinks, 'Vartype','unequal', 'tail', 'right');

% make a figure
f=figure('color','white');
set(gcf, 'Position', [10 10 600 600]);
myColors = GiveMeColors('RFPU');
linkNames = {'rich', 'feeder', 'peripheral', 'unconnected'};

for j=1:length(dataCell)
    J.(linkNames{j}) = dataCell{j};
end

violins = violinplot(J);
set(gcf,'color','w');
xtickangle(30);
set(gca,'fontsize', 18);

if isCorrected
    ylabel({'Spatially-corrected', 'correlated gene expression'})
else
    ylabel({'Correlated gene expression'})
end
ylim([-1 1])


for i=1:size(myColors,1)
    violins(1,i).ViolinColor = myColors(i,:);
    violins(1,i).EdgeColor = myColors(i,:);
    violins(1,i).BoxColor = [.25 .25 .25];
    
end

%t = stats.tstat;
end
