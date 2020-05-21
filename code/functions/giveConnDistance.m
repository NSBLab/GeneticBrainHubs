function distMatr = giveConnDistance(parcellation, tractography, groupAdjlog)

if strcmp(parcellation, 'HCP')
    conFile = 'HCPMMP1ANDfslatlas20';
    cortex = [1:180,191:370];
elseif strcmp(parcellation, 'random500')
    conFile = 'random500ANDfslatlas20';
    cortex = [1:250,261:510];
elseif strcmp(parcellation, 'cust250')
    conFile = 'custom500ANDfslatlas20';
    cortex = [1:250,261:510];
end

load(sprintf('HCP_%s_%s_NOSIFT_length_structnets.mat', conFile, tractography));
% take the average ofconnection lengths to be used for distance effect
% checking.
conLength = zeros(size(ADJS{1},1),size(ADJS{1},1),length(ADJS));
for sub=1:length(ADJS)
    conLength(:,:,sub) = ADJS{sub};
end
conLength(:,:,299) = [];
avLength = mean(conLength,3); % get the mean of lengths
distMatr = avLength(cortex,cortex); % select cortical regions
distMatr = distMatr.*groupAdjlog; % get values only for existing links
distMatr(groupAdjlog==0) = NaN;
 
end

