function run_extract_null_eQTLs(numNull)

cd /projects/kg98/aurina/GeneticBrainHubs
% load paths
%-------------------------------------------------------------------------------
fprintf(1,'Adding all subdirectories to the Matlab path...');
% Add paths required for the project (ignoring hidden, including version control)
files = dir;
directories = files([files.isdir]);
directories(strmatch('.',{files([files.isdir]).name})) = []; % remove hidden
if isfield(directories,'folder')
    paths = arrayfun(@(x)fullfile(directories(x).folder,directories(x).name),1:length(directories),'UniformOutput',false);
else
    paths = arrayfun(@(x)fullfile(pwd,directories(x).name),1:length(directories),'UniformOutput',false);
end
for j = 1:length(paths)
    addpath(genpath(paths{j}))
end
fprintf(1,' Added.\n');
%-------------------------------------------------------------------------------

% run for HCP
[HCP_pORA_eQTL, HCP_genes_ENTREZ] = extract_null_eQTLs('HCP', 'POS');
fileNameHCP = sprintf('/projects/kg98/aurina/GeneticBrainHubs/data/reeqtls/permutations/HCP_POS_null_%d.mat', numNull); 
save(fileNameHCP,'HCP_pORA_eQTL', 'HCP_genes_ENTREZ');

% run for GenCog
[GenCog_pORA_eQTL, GenCog_genes_ENTREZ] = extract_null_eQTLs('GenCog', 'POS');
fileNameGenCog = sprintf('/projects/kg98/aurina/GeneticBrainHubs/data/reeqtls/permutations/GenCog_POS_null_%d.mat', numNull); 
save(fileNameGenCog,'GenCog_pORA_eQTL', 'GenCog_genes_ENTREZ');

end