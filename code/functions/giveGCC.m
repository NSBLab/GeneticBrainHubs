% Gene scorring
% for every gene there is a region x region matrix defining the
% contribution of that particular gene

function GCC = giveGCC(parcelExpression, selectGenes, FitCurve)

expMatrix = parcelExpression(:,selectGenes+1);
numRegions = size(parcelExpression,1);
numGenes = size(expMatrix,2);
% calculate gene contribution score for each gene in every region pair using distance corrected data
GCC = zeros(numRegions, numRegions, numGenes);

for i=1:numRegions
    zscoredVal1 = zscore(expMatrix(i,:));
    
    for j=i+1:numRegions
        zscoredVal2 = zscore(expMatrix(j,:));
        
        for g=1:numGenes
            
            GCC(i,j,g) = zscoredVal1(g).*zscoredVal2(g) - FitCurve(i,j);
            
        end
        
    end
    
end

end


