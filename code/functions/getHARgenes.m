% import HAR genes from Doan, 2016. https://doi.org/10.1016/j.cell.2016.08.071, select unique genes
HARunique = (importHARgenes('HAR_genes.xlsx')); 
% save the list of genes for converting names to entrez IDs
writetable(HARunique, 'data/HARgenes/HAR_uniqueGenes.txt')
