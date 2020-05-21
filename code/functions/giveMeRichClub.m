% Rich club in the whole brain

function [G,  mdist, PhiNormMean, dMiddle] = giveMeRichClub(matrices, COG, groupMatrixType, densThreshold, giveRC, CVmeasure, consThr, numIter, numRepeats, whatTypeNetwork ,whatNullModel,yVals)

if nargin < 3
    groupMatrixType = 'consistency';
    densThreshold = 0.8; 
    giveRC = false; 
    fprintf('Making consistency based group matrix BY DEFAULT\n')
end

if nargin < 10
    if giveRC
        numIter = 50;
        numRepeats = 100;
        %     whatTypeNetwork = 'wu';
        %     whatNullModel = 'geometry';
        whatTypeNetwork = 'bu';
        whatNullModel = 'randmio_und';
        yVals = [0.95 1.5];
        if giveRC
            fprintf('Calculating RC using %d iterations with %d repeats on %s network with %s null model BY DEFAULT\n', numIter, numRepeats, whatTypeNetwork, whatNullModel)
        end
    end
end

if nargin < 12
    yVals = [0.95 1.5];
end



numNodes = size(matrices{1},1);
numSubjects = length(matrices);

A = zeros(numNodes, numNodes,numSubjects);
dist = zeros(numNodes, numNodes,numSubjects);
for m=1:numSubjects
    A(:,:,m) = matrices{m};
    coords = COG{m};
    dist(:,:,m) = pdist2(coords,coords);
end
mCon = mean(A,3);
mdist = mean(dist,3);

% get the weights for only the existing edges
nanA = A;
nanA(nanA==0) = NaN;
avWeight = nanmean(nanA,3);
avWeight(isnan(avWeight)) = 0;

% make a group matrix
if strcmp(groupMatrixType, 'lengthCV')
    if numNodes==180
        hemiid = ones(numNodes,1);
    else
        
    hemiid = zeros(numNodes,1);
    hemiid(1:numNodes/2) = 1;
    hemiid(numNodes/2+1: numNodes) = 2;
    end
    
    Gr = fcn_group_average(A,mdist,hemiid);
    G = Gr.*avWeight;
    
elseif strcmp(groupMatrixType, 'CVmeasure')
    Gr = giveMeGroupAdj_variance_AA(matrices, densThreshold, CVmeasure, consThr);
    Gr = logical(Gr); 
    G = Gr.*avWeight;

elseif strcmp(groupMatrixType, 'consistency')
    G = giveMeGroupAdj_consistency(matrices, densThreshold);
end

if giveRC
    if strcmp(whatNullModel, 'shuffleWeights') && max(G(:))>100
        G = log(G); 
        G(isinf(G)) = 0; 
    end
[~, dMiddle, PhiNormMean] = PlotRichClub(G,mdist,whatTypeNetwork,whatNullModel,numIter,numRepeats,yVals);
end 

end