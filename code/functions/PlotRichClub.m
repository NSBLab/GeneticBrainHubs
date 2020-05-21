function [dMiddleNorm, dMiddle, PhiNormMean, PhiTrue, PhiRand, isSig]=PlotRichClub(Adj,DISTmean, WhatTypeNetwork,whatNullModel,numIter,numRepeats,yVals, whatDistribution, colorOut, colorIn)
%-------------------------------------------------------------------------------
% Loads in data from a computation using networkRC,
% Plots RC curve and degree distribution
%-------------------------------------------------------------------------------
%colorInbarCo = [198 219 239]/255; 
colorInbarCo = colorOut; 
doRankSum = false; % Welch's t-test or ranksum test?
meanOrMedian = 'mean'; % mean distance

%-------------------------------------------------------------------------------
% Load in the data
%-------------------------------------------------------------------------------


if strcmp(WhatTypeNetwork, 'bu') || strcmp(WhatTypeNetwork, 'bd')
    Adj = logical(Adj);
elseif strcmp(WhatTypeNetwork, 'wu') || strcmp(WhatTypeNetwork, 'wd')
%     Adj = log(Adj);
%     Adj(isinf(Adj)) = 0;
end

if strcmp(WhatTypeNetwork, 'bd')
    [~,~,deg] = degrees_dir(Adj);
else
    deg = degrees_und(Adj);
end

if strcmp(WhatTypeNetwork, 'wuStrength') || strcmp(WhatTypeNetwork, 'wuStrengthBins')
    kmax = 100;
else
    kmax = max(deg);
end

pThreshold = 0.05;
[PhiNorm,PhiTrue,PhiRand] = RichClubPhiNorm(Adj,kmax, numIter,numRepeats,WhatTypeNetwork,whatNullModel, DISTmean); %, doBins);
pValues = zeros(kmax,1);


for i = 1:kmax
    pValues(i) = mean(PhiTrue(i) <= PhiRand(:,i));
end

% Significant phi
isSig = (pValues <= pThreshold);
PhiNormMean = zeros(size(PhiTrue));

for i = 1:length(PhiTrue)
    PhiNormMean(i) = PhiTrue(i)/mean(PhiRand(:,i));
end

% %===============================================================================
% % Distance curve
% %===============================================================================
kRange = min(deg):max(deg);
% Add a baseline (1)
dAll = cell(length(kRange),2);
pVals = zeros(length(kRange),1);
for i = 1:length(kRange)
    isHub = double(deg > kRange(i));
    isRich = isHub'*isHub;
    isRich(eye(logical(size(isRich)))) = 0;
    notRich = ~isRich;
    notRich(eye(logical(size(isRich)))) = 0;
    dAll{i,1} = DISTmean(Adj & isRich);
    dAll{i,2} = DISTmean(Adj & notRich);
    
    % Difference
    if all(isnan(dAll{i,1})) || all(isnan(dAll{i,2}))
        pVals(i) = NaN;
    else
        if doRankSum
            pVals(i) = ranksum(dAll{i,1},dAll{i,2},'tail','right');
        else
            [~,pVals(i)] = ttest2(dAll{i,1},dAll{i,2},'Tail','right','Vartype','unequal');
        end
    end
end

switch meanOrMedian
    case 'mean'
        dMiddle = cellfun(@nanmean,dAll);
    case 'median'
        dMiddle = cellfun(@nanmedian,dAll);
end
dMiddleNorm = dMiddle(:,1);
for j=1
    
    figure('color','w');
    axMain = subplot(5,3,7:15);hold on;
    axMain.Position = [0.13 0.18 0.775 0.377];
    plot([min(deg),max(deg)],ones(2,1),':','color','k','LineWidth',3)
    %===============================================================================
    % Add rich-club curve:
    %===============================================================================
    plot(PhiNormMean, '-','Color', colorOut ,'LineWidth',2.5);
    % Significance at p = 0.05 as circles
    plot(find(isSig),PhiNormMean(isSig),'o','MarkerEdgeColor',colorOut,...
        'MarkerFaceColor',colorIn,'LineWidth',1.5,'MarkerSize',6);
    % Add shaded rectangle
    %h_rect = rectangle('Position',[30,ylimNow(1),60,ylimNow(2)-ylimNow(1)+.1],...
    % 'EdgeColor','none','FaceColor',ones(3,1)*0.90);
    %uistack(h_rect,'bottom');
    
    % Set axis limits, legend, labels:
    yVals = [nanmin(PhiNormMean(isSig))-0.5*nanstd(PhiNormMean) nanmax(PhiNormMean(isSig))+0.5*nanstd(PhiNormMean)]; 
    axMain.YLim = yVals;
    xlabel('Node degree, k','fontsize',18);
    ylabel('\Phi_{norm}','fontsize',18);
    get(gca, 'XTick');
    set(gca, 'FontSize', 18)
    box off;
    %legend('boxoff')
    box('off')
    axMain.XLim = [min(deg)-0.5,max(deg)+2];
    
    
    axDD_type = subplot(5,3,1:6);
    switch whatDistribution
        
        case 'barCount'
            
            N = arrayfun(@(x)sum(deg==x),kRange);
            % this plots for every degree 
            %bar(kRange,N,'EdgeColor',colorOut,'FaceColor',colorInbarCo);
           
            % this version averages bin counts into ranges (like historgam)
            
            numBins = 6; % 2 works as well
            u=1;
            w=1;
            while u<=length(N)
                Nbin(w) = sum(N(u:u+numBins));
                kRangebin(w) = mean(kRange(u:u+numBins));
                u=u+numBins+1;
                w=w+1;
            end
            
            figure;
            bar(kRangebin,Nbin,'EdgeColor',colorOut,'FaceColor',colorInbarCo);
            
        case 'histogram'
            
            histogram(deg,35,'EdgeColor',[1 1 1],'FaceColor',colorInbarCo, 'FaceAlpha',1); 
            
        case 'kernel'
            
            pdSix = fitdist(deg','Kernel','BandWidth',5);
            x = 0:1:max(deg)+5;
            ySix = pdf(pdSix,x);
            plot_shaded(x,ySix, 'Color', colorOut)
            %plot(x,ySix,'k-','LineWidth',2)
    end
    
    xticks([]); box off;
    xlim([min(deg)-0.5 max(deg)+2]);
    %ylim([0 max(a.YData)+1]);
    axDD_type.Position = [0.13 0.6 0.775 0.3];
    
    ylabel('Frequency','fontsize',18);
    get(gca, 'YTick');
    set(gca, 'FontSize', 18)
    %axDD_type.Position = [0.1300    0.776    0.775    0.149];
    set(gcf, 'Position', [500 500 750 550])
    %pos=get(axMain,'Position');
    
    box off;
    %     if j==1
    %     figureName = sprintf('results/figures/RC_%s_%dnodes_dens%d.png', WhatTypeNetwork, size(Adj,1), density_und(Adj)*100);
    %     else
    %     figureName = sprintf('results/figures/DIST_%s_%dnodes_dens%d.png', WhatTypeNetwork, size(Adj,1), density_und(Adj)*100);
    %     end
    %     print(gcf,figureName,'-dpng','-r600');
end
end


