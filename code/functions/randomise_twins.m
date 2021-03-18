% This script is randomising familial relationships (in the same order across edges)
% while keeping the pairings of twins the same by swapping pairs;
% MZ/DZ twins are randomised together;
 
% The randomised connectivity and covariate values are saved as a separate file 
% and will be used to run heritability on each file separately

% -------------
% LOAD data
% -------------
load('twinEdges_HCP_iFOD2_FA_strength20.mat')
load('twinCovariatesDWI.mat')

% -------------
% select only MZ and DZ twins, exclude siblings and merge MZ and DZ twins together
% -------------
Output_MZ = Output_MZ(:,1:2,:);
Output_DZ = Output_DZ(:,1:2,:);

MZ_ID = MZ_ID(:,1:2); 
DZ_ID = DZ_ID(:,1:2);
ID_ALL = vertcat(MZ_ID, DZ_ID);

MZ_age = MZ_age(:,1:2); 
DZ_age = DZ_age(:,1:2);
AGE_ALL = vertcat(MZ_age, DZ_age);

MZ_sex = MZ_sex(:,1:2); 
DZ_sex = DZ_sex(:,1:2);
SEX_ALL = vertcat(MZ_sex, DZ_sex);

% save only twins for a separate run of the analyses
save('data/heritability/twinEdges_HCP_iFOD2_FA_strength20_only_twin.mat', 'Output_MZ', 'Output_DZ', 'groupAdjlog')
save('data/heritability/twinCovariatesDWI_only_twin.mat', 'MZ_ID', 'DZ_ID','MZ_age', 'DZ_age', 'MZ_sex', 'DZ_sex')

Output_ALL = vertcat(Output_MZ, Output_DZ);

% randomise 20 times and save as separate files
for i=1:20
    % randomise the order across subjects
    idx = randperm(size(Output_ALL,1));
    
    % -------------
    % CONNECTIVITY flip connectivity values for all edges
    % -------------
    Output_ALLnew = Output_ALL(:,1,:);
    Output_ALLnew(idx,:,:) = Output_ALLnew(flip(idx),:,:);
    
    % replace old values with new
    Output_ALL_rand = Output_ALL;
    Output_ALL_rand(:,1,:) = Output_ALLnew;
    
    % -------------
    % ID flip values for all edges
    % -------------
    ID_ALLnew = ID_ALL(:,1,:);
    ID_ALLnew(idx,:,:) = ID_ALLnew(flip(idx),:,:);
    
    % replace old values with new
    ID_ALLnew_rand = ID_ALL;
    ID_ALLnew_rand(:,1,:) = ID_ALLnew;
    
    % -------------
    % AGE flip connectivity values for all edges
    % -------------
    AGE_ALLnew = AGE_ALL(:,1,:);
    AGE_ALLnew(idx,:,:) = AGE_ALLnew(flip(idx),:,:);
    
    % replace old values with new
    AGE_ALLnew_rand = AGE_ALL;
    AGE_ALLnew_rand(:,1,:) = AGE_ALLnew;
    
    % -------------
    % SEX flip connectivity values for all edges
    % -------------
    SEX_ALLnew = SEX_ALL(:,1,:);
    SEX_ALLnew(idx,:,:) = SEX_ALLnew(flip(idx),:,:);
    
    % replace old values with new
    SEX_ALLnew_rand = SEX_ALL;
    SEX_ALLnew_rand(:,1,:) = SEX_ALLnew;
    
    
    % separate subjects into MZ and DZ twins at random
    Output_MZ_rand = Output_ALL_rand(1:size(Output_MZ,1),:,:);
    Output_DZ_rand = Output_ALL_rand(size(Output_MZ,1)+1:end,:,:);
    
    % same with covariates
    MZ_ID_rand = ID_ALLnew_rand(1:size(Output_MZ,1),:,:);
    DZ_ID_rand = ID_ALLnew_rand(size(Output_MZ,1)+1:end,:,:);
    
    MZ_age_rand = AGE_ALLnew_rand(1:size(Output_MZ,1),:,:);
    DZ_age_rand = AGE_ALLnew_rand(size(Output_MZ,1)+1:end,:,:);
    
    MZ_sex_rand = SEX_ALLnew_rand(1:size(Output_MZ,1),:,:);
    DZ_sex_rand = SEX_ALLnew_rand(size(Output_MZ,1)+1:end,:,:);
    
    filename_con = sprintf('data/heritability/twinEdges_HCP_iFOD2_FA_strength20_rand%d.mat', i);
    save(filename_con, 'Output_MZ_rand', 'Output_DZ_rand');
    
    filename_cov = sprintf('data/heritability/twinCovariatesDWI_rand%d.mat', i);
    save(filename_cov, 'MZ_ID_rand', 'DZ_ID_rand', 'MZ_age_rand', 'DZ_age_rand','MZ_sex_rand', 'DZ_sex_rand');
    
end