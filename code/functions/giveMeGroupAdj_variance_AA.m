function [groupAdj, consist] = giveMeGroupAdj_variance_AA(connectomes, dens, CVmeasure, consThr)

M = zeros(size(connectomes{1},1),size(connectomes{1},2),size(connectomes,2));
nSubs = size(connectomes,2);
d = zeros(nSubs,1);
for i=1:nSubs
    M(:,:,i) = connectomes{i};
    d(i) = density_und(connectomes{i});
end
dMean = mean(d);
% if density is specified, use that value
if nargin >1
    dMean = dens;
end
    [groupAdj, consist] = threshold_consistency_AA(M, dMean, CVmeasure, consThr);

end
