function [rAdj] = shuffleWeights_herit(Adj)
        % Ben Fulcher, 2014-12-01
        % Shuffles weights, keeping the topology fixed
        
        % 1. Find
        % Get all elements of link data where a connection exists:
        allActualLinks = Adj(~isnan(Adj));
        
        % Shuffle them:
        allActualLinksDataShuffled = allActualLinks(randperm(length(allActualLinks)));
        
        % Put them back in the matrix
        rAdj = ~isnan(Adj)+0;
        rAdj(rAdj==1) = allActualLinksDataShuffled;
        
        % Not relevant to this method:
end