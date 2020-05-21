function [BBmean,microDist,distNorm, mpc_bb, mp, BBskew] = getBBdata(parc, yregress)

% get regional proviles, mean values, skewness distances in the gradient space.

%load(sprintf('bb_%s.mat', parc))

% load corresponding connectivity matrix
if strcmp(parc, 'HCP')
    %fileAdj = 'HCP_360_20.mat';
    fileMicro = 'glasser';
elseif strcmp(parc, 'custom200')
    % fileAdj = 'HCP_200_25.mat';
    fileMicro = 'custom200';
elseif strcmp(parc, 'custom500')
    %fileAdj = 'HCP_500_15.mat';
    fileMicro = 'custom500';
elseif strcmp(parc, 'random500')
    %fileAdj = 'HCP_500_15.mat';
    fileMicro = 'random500';
elseif strcmp(parc, 'random200')
    %fileAdj = 'HCP_500_15.mat';
    fileMicro = 'random200';
end

%load(fileAdj)

if yregress
    fileProfiles = sprintf('BB_profiles_Yregress_%s.csv', fileMicro); 
else
    fileProfiles = sprintf('BB_profiles_%s.csv', fileMicro);     
end

mp = importBBprofiles(fileProfiles)'; 
BBlabels = importBBreg(fileProfiles);  
    
% reorder regions for HCP parcellation
if strcmp(parc, 'HCP') || strcmp(parc, 'HCPVcorr')
    
    EXPlabels = importExpRegLabels('newlabels_HCPMMP1ANDfslatlas20.txt');
    EXPlabels = EXPlabels(:,2);
    [~, ~, indBB] = intersect(EXPlabels,BBlabels, 'stable');
    % reorder regions
    mp = mp(:,indBB); 
end
    
    % mean across layers
    BBmean = mean(mp,1);
    BBskew = skewness(mp); 
   

    microDist = squareform(pdist(BBmean'));
    distNorm = microDist./(max(microDist(:)));
    
    % create MPC matrix from partial correlation
    R = partialcorr(mp, mean(mp,2));
    % remove negative values
    R(R<0) = 0;
    % log transformation
    mpc_bb = 0.5 * log( (1+R) ./ (1-R));
    %mpc_bb(isnan(mpc_bb)) = 0; 
    mpc_bb(isinf(mpc_bb)) = 0;
    mpc_bb(mpc_bb==0) = NaN; 
    

end