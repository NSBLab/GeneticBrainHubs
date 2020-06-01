function plot_modellingCDF(Adj, Adj_mod, Adj_dist, numNetworks)

NetProp{1} = sum(Adj);
NetProp{2} = betweenness_bin(Adj);
NetProp{3} = clustering_coef_bu(Adj);
NetProp{4} = Adj_dist(triu(Adj,1) > 0);

SimProp = cell(4,numNetworks);

for j = 1:numNetworks
    
    SimProp{1,j} = sum(Adj_mod{j});
    SimProp{2,j} = betweenness_bin(Adj_mod{j});
    SimProp{3,j} = clustering_coef_bu(Adj_mod{j});
    SimProp{4,j} = Adj_dist(triu(Adj_mod{j},1) > 0);
    
end

numMeasures = size(NetProp,2);

for i = 1:numMeasures
    
    subplot(1,numMeasures,i)
    
    [f,x] = ecdf(NetProp{i});
    plot(x,f,'Color','k')
    hold on
    for j = 1:numNetworks
        [f,x] = ecdf(SimProp{i,j});
        plot(x,f,'Color','r')
    end
    
end


end
