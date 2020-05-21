function rAdj = randomiseWeights(herit, groupAdjlog)
% Ben Fulcher, 2014-12-01
% Shuffles weights, keeping the topology fixed

% 1. Find
% Get all elements of link data where a connection exists:
allActualLinks = herit(groupAdjlog==1);

% Shuffle them:
allActualLinksDataShuffled = allActualLinks(randperm(length(allActualLinks)));

% Put them back in the matrix
rAdj = logical(groupAdjlog)+0;
rAdj(groupAdjlog==1) = allActualLinksDataShuffled;
rAdj(groupAdjlog==0) = NaN; 
end