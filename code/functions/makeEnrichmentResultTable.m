
function makeEnrichmentResultTable()
% Compile ermineJ results and keep only signifficant values after FDR correction

% this is output file from ermineJ
type = 'Rich_LRcortex_2.000000e-01perc_3.000000e-01_strength_k105_HCP_RNAseqprobes_NWDC_pearson_03_05_2021'; 
fileINname = sprintf('%s.txt', type);  

ermineJResultsBPcon = ReadInErmineJ(fileINname); 
CONall = vertcat(ermineJResultsBPcon); 
CONall.corr_pval = mafdr(CONall.pval,'BHFDR',true);
CONall = CONall(1:100,:);
GOtableCON = table(); 
GOtableCON.GOcategory = CONall.GOID;
GOtableCON.Description = CONall.GOName;
GOtableCON.NumGenes = CONall.numGenes; 
num_dig = 15;
GOtableCON.Pval = round(CONall.pval*(10^num_dig))/(10^num_dig);
GOtableCON.Pval_corr = round(CONall.corr_pval*(10^num_dig))/(10^num_dig);

% save file as .csv
cd('data/enrichment'); 
fileOUTname = sprintf('ermineJresults%s.csv', type); 
writetable(GOtableCON,fileOUTname)

end




