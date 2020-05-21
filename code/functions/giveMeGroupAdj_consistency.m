
function [groupAdj, mconsist] = giveMeGroupAdj_consistency(connectomes, threshold)

% make masks for within and between hemisheres
% make 3D matrix and calculate mean distances
numSubjects = length(connectomes); 
numNodes = size(connectomes{1},1); 

%distance = zeros(numNodes,numNodes, numSubjects); 
adj = zeros(numNodes,numNodes, numSubjects); 
consist = zeros(numNodes,numNodes, numSubjects);  

for s = 1:numSubjects
    %distance(:,:,s) = distances{s}; 
    adj(:,:,s) = connectomes{s}; 
    consist(:,:,s) = logical(connectomes{s}); 
end

%%distance(distance==0) = NaN; 
%mLength = nanmean(distance,3); 

% replace zeros with NaN in order to take the average of only existing
% values (and not average non-existing zeros)
adj(adj==0) = NaN; 
madj = nanmean(adj,3); 

mconsist = mean(consist,3); 

Adj = madj;
mask = mconsist>=threshold; 
groupAdj = Adj.*mask; 
groupAdj(isnan(groupAdj)) = 0; 

%groupDist = mLength.*mask; 
%groupDist(isnan(groupDist)) = 0; 

end





