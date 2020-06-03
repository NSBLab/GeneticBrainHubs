function plot_modellingCDF(Adj, Adj_mod, Adj_dist, numNetworks)

[rgb_colorMatrix] = GiveMeColors('4orangepurple'); 
% make a little lighter for model lines
rgb_colorMatrixLines = brighten(rgb_colorMatrix, 0.3); 
rgb_colorMatrixReal = brighten(rgb_colorMatrix, -0.3); 

NetProp{1} = sum(Adj);
measureNames{1} = {'Node degree, k', 'k'}; 
NetProp{2} = betweenness_bin(Adj);
measureNames{2} = {'Betweenness, b' ,'b'}; 
NetProp{3} = clustering_coef_bu(Adj);
measureNames{3} = {'Clustering, c', 'c'}; 
NetProp{4} = Adj_dist(triu(Adj,1) > 0);
measureNames{4} = {'Edge length, d', 'd'}; 

SimProp = cell(4,numNetworks);

for j = 1:numNetworks
    
    SimProp{1,j} = sum(Adj_mod{j});
    SimProp{2,j} = betweenness_bin(Adj_mod{j});
    SimProp{3,j} = clustering_coef_bu(Adj_mod{j});
    SimProp{4,j} = Adj_dist(triu(Adj_mod{j},1) > 0);
    
end

numMeasures = size(NetProp,2);

si = [33 37 41 45]; 
for i = 1:numMeasures
    
    subplot(12,numMeasures,si+i-1); 
    for j = 1:numNetworks
        [f,x] = ecdf(SimProp{i,j});
        plot(x,f,'Color',rgb_colorMatrixLines(i,:))
        axis square
        box off
        hold on
    end
    % plot measures for real data on top of models    
    [f,x] = ecdf(NetProp{i});
    plot(x,f,'Color',rgb_colorMatrixReal(i,:),'LineWidth', 4);
    xlabel(measureNames{i}{1})
    ylabel(sprintf('F(%s)', measureNames{i}{2}))
    set(gca,'FontSize',16)
    xlim([0 max(NetProp{i})])
    if strcmp(measureNames{i}{1}, 'Betweenness, b')
    xticks([0 500 1000])
    end
end


end
