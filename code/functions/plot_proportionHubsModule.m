function F = plot_proportionHubsModule(netassignments, nodeData,plotOptions)

colorOut = plotOptions.colorOut;         
modulesColor = GiveMeColors('funcModules'); 

k=linspace(1,max(nodeData), max(nodeData)); 
numMod = max(unique(netassignments)); 

Hp = zeros(numMod, max(k)-1); 
for i=1:max(k)-1
    isHub = nodeData>k(i); 
    numHubs = length(find(isHub)); 
    %proportion of hubs in each module
    for m=1:numMod
        isMod = netassignments==m; 
        isModHub = isMod&isHub'; 
        numModHub = length(find(isModHub)); 
        Hp(m,i) = numModHub./numHubs; 
    end
end

% plot degree distribution
F = figure('color','w');
set(gcf, 'Position', [500 500 750 550])
sp=subplot(5,3,1:6);
histogram(nodeData,35,'EdgeColor',[1 1 1],'FaceColor',colorOut, 'FaceAlpha',1);
xlim([min(nodeData)-0.8,max(nodeData)-0.5]);
xticks([]); box off;
ylabel('Frequency')
%xlabel('Node degree, k')
get(gca, 'YTick');
set(gca, 'FontSize', 20)
hold on;

% plot the proportion
sp=subplot(5,3,7:15);
H = bar(Hp','stacked','BarWidth', 1);
ylim([0 1]); box off; 
xlim([min(nodeData)-0.8,max(nodeData)-0.5]);
ylabel({'Proportion of nodes', 'with degree>k'})
xlabel('Node degree, k')
set(gca, 'FontSize', 20)

% change color for each module
for i = 1:size(modulesColor,1)
    H(i).FaceColor = 'flat';
    H(i).CData = modulesColor(i,:); 
end



end

