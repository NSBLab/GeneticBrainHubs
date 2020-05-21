function ComNorm = normCommunicability(Adj, whatType, numIter, numRepeats, whatNullModel)

if strcmp(whatType, 'bin')
    weighted = 0;
elseif strcmp(whatType, 'wei')
    weighted = 1;
end

ComTrue = communicability(Adj,weighted);

switch whatNullModel
    case 'randmio_und'
        f_rand_null = @randmio_und;
    case 'shuffleWeights'
        f_rand_null = @f_shuffleWeights;
end

ComRand = zeros(numRepeats, size(Adj,1), size(Adj,1));

for i=1:numRepeats
    
    timer = tic;
    fprintf(1,'[%u/%u] Rewiring each link %u times...\n',i,numRepeats,numIter);
    
    [Adj_rand] = f_rand_null(Adj,numIter); % Random graph with preserved in/out degree distribution
    ComRand(i,:,:) = communicability(Adj_rand,weighted);
    
    if i==1 || mod(i,numRepeats/10)==0
        fprintf(1,'Approx. %s remaining...\n',BF_thetime(toc(timer)/i*(numRepeats-i)));
    end
    
end

% calculate the mean of non-nan communicability values in random networks
isNaN = isnan(ComRand); 
isINF = isinf(ComRand);

% label NaN and inf values as NaNs
badLinks = isNaN | isINF; 
ComRandKEEP = ComRand;
ComRandKEEP(badLinks==1) = NaN; 

% get nanmean ignoring those values
ComRandMean = squeeze(nanmean(ComRandKEEP,1)); 
% divide empirical by random and get normalised communicability
ComNorm = ComTrue./ComRandMean; 
% ------------------------------------------------------------------------------
% Extra functions:
    function [rAdj,numRewirings] = f_shuffleWeights(Adj,numIter);
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
        numRewirings = 0;
    end
% ------------------------------------------------------------------------------
end
