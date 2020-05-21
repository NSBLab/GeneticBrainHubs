function [RvsF, FvsP, dataCellALL, xThresholds,f] = plot_distanceViolin(pairwiseMeasure, distMatr, groupAdjlog, nodeData, khub, numThr, yTtile)

% select one half and NaN values for non-ezisting edges just in case it was
% not done before
pairwiseMeasure = maskuHalf(pairwiseMeasure);
pairwiseMeasure(groupAdjlog==0) = NaN;


[xThresholds] = BF_PlotQuantiles(distMatr(:), pairwiseMeasure(:),numThr,0,0);
[Y] = discretize(distMatr,xThresholds);

f=figure('color','white');
set(gcf, 'Position', [10 10 1700 600]);
%set(gcf, 'Position', [10 10 1200 500]);
myColors = GiveMeColors('RFPU');
linkNames = {'rich', 'feeder', 'peripheral'};

for b=1:max(Y(:))
    
    maskBin = Y==b; % all links in this distance bin
    AdjCon = groupAdjlog.*Y==b; % connected links in this distance bin
    
    [dataCell,pRF,statsRF,pFP,statsFP] = AA_HERIT_RFPU(pairwiseMeasure, AdjCon, nodeData, khub);
    
    RvsF(b,2) = pRF;
    RvsF(b,1) = statsRF.tstat;
    
    FvsP(b,2) = pFP;
    FvsP(b,1) = statsFP.tstat;
    % use data from dataCell for making distributions
    for j=1:length(dataCell)
        J.(linkNames{j}) = dataCell{j};
    end
    
    dataCellALL{b} = dataCell;
    
    subplot(1,numThr-1,b);
    set(gca,'FontSize',18)
    violins = violinplot(J);
    set(gcf,'color','w');
    xtickangle(30);
    
    switch yTtile
        case 'Heritability'
        ylim([0 1])
        yLAB = 'Heritability';
        case 'CGE'
        yLAB = {'Spatially-corrected', 'correlated gene expression'};
        ylim([-0.8 0.8])   
    end
    ylabel(yLAB)
    
    for i=1:size(myColors,1)-1
        violins(1,i).ViolinColor = myColors(i,:);
        violins(1,i).EdgeColor = myColors(i,:);
        violins(1,i).BoxColor = [.25 .25 .25];
        %violins(1,i).ShowMean = 1;
    end
    
   
    for l=1:size(dataCell,1)
        
        A = strcat(linkNames{l}, {' '}, {'('}, num2str(length(dataCell{l})),{')'});
        violinLabels{l} = A{1};
        
    end
    
    xticklabels(violinLabels);
    xtickangle(30)
    hold on;
    
end