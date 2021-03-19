function plot_onlyHist(nodeData, plotOptions)

whatDistribution = plotOptions.whatDistribution; 
colorOut = plotOptions.colorOut; 
colorIn = plotOptions.colIn; 

sortK = sort(nodeData,'descend');
maxK = sortK(2); % Up to the second-highest k

% All k are a bin
kr = min(nodeData):maxK;

krAll = min(nodeData):max(nodeData);


figure('color','w');
sp=subplot(5,3,1:6);

switch whatDistribution
    case 'barCount'
        
        N = arrayfun(@(x)sum(nodeData==x),krAll);
        bar(krAll,N,'EdgeColor',colorOut,'FaceColor',colorOut);
        
    case 'histogram'
        
        histogram(nodeData,35,'EdgeColor',[1 1 1],'FaceColor',colorOut, 'FaceAlpha',1);
        
    case 'kernel'
        
        pdSix = fitdist(nodeData','Kernel','BandWidth',3);
        x = 0:1:max(nodeData)+5;
        ySix = pdf(pdSix,x);
        plot_shaded(x,ySix, 'Color', colorOut)
        %plot(x,ySix,'k-','LineWidth',2)
end

xlim([min(nodeData)-0.8,max(nodeData)+0.8]);
xticks([]); box off;
ylabel('Frequency', 'FontSize', 18)
get(gca, 'YTick');
set(gca, 'FontSize', 18)
sp=subplot(5,3,7:15);
hold on;

end

