for i=1:size(coordinates,2)
    A(i,:,:,:) = coordinates{i}; 
end
% mean over subjects
A = A(:,1:180,:); 
mA = squeeze(mean(A, 1)); 

% for each coorditate, get max difference in values
for k=1:3
    for j=1:size(mA,1)
        for l=j+1:size(mA,1)
            N(j,l) = abs(mA(j,k) - mA(l,k)); 
        end
    end
    maxK(k) = max(N(:)); 
end
