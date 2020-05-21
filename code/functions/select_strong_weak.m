function mask = select_strong_weak(dataMatrix, n)

%---INPUTS:
% dataMatrix, the input data matrix
% n, the number of links to select in each group

%---OUTPUT:
% mask - matrix, indicating TOP x number of strongest and weakest links in the matrix

ind = zeros(length(find(~isnan(maskuHalf(dataMatrix)))),2); 
ind(:,1) = find(~isnan(maskuHalf(dataMatrix))); 
ind(:,2) = dataMatrix(~isnan(maskuHalf(dataMatrix(:)))); 

indS = sortrows(ind,-2); 

% select TOP links
indTOP = indS(1:n,:); 

% select BOTTOM links
% select minimum absolute values (edges showing weakest relationships)
[~, I] = min(abs(indS(:,2))); 
indBOT = indS(I-n/2:I+n/2-1,:); 

% % put vallues back to matrix
mask = zeros(size(dataMatrix(:))); 
mask(indTOP(:,1)) = 2; 
mask(indBOT(:,1)) = 1; 
mask = reshape(mask,size(dataMatrix)); 



end