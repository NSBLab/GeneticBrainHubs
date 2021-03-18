function [rAdj] = shuffleWeights(Adj)
        % Ben Fulcher, 2014-12-01
        % Shuffles weights, keeping the topology fixed
        
        % 1. Find
        % Get all elements of link data where a connection exists:
        allActualLinks = Adj(Adj~=0);
        
        % Shuffle them:
        allActualLinksDataShuffled = allActualLinks(randperm(length(allActualLinks)));
        
        % Put them back in the matrix
        rAdj = logical(Adj)+0;
        rAdj(Adj~=0) = allActualLinksDataShuffled;
        
        % Not relevant to this method:
end