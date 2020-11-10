%--------------------------------------------------% 
% This script runs all the functions to reproduce figures in the manuscript
%--------------------------------------------------%

setupPaths()
%--------------------------------------------------%
% Figure 2 - main text: heritability
%--------------------------------------------------%
Figure2()

%--------------------------------------------------%
% Figure 3 - main text: correlated gene expression
%--------------------------------------------------%
Figure3()

%--------------------------------------------------%
% Figure 4 - main text: generative modelling
%--------------------------------------------------%
Figure4('energy')

%--------------------------------------------------%
% Figure S1 - supplementary: topological properties 
% of the connectome
%--------------------------------------------------%
FigureS1()

%--------------------------------------------------%
% Figure S2-S5 - supplementary: heritability-related analyses
%--------------------------------------------------%
FigureS2_S5()

%--------------------------------------------------%
% Figure S6-S11 - supplementary: correlated gene 
% expression related analyses
%--------------------------------------------------%
FigureS6_S11()

%--------------------------------------------------%
% Figure S12 - supplementary: microstructure related analysis
%--------------------------------------------------%
FigureS12()

%--------------------------------------------------%
% Figure S13 - supplementary: generative modelling 
%--------------------------------------------------%
Figure4('degCorr')

%--------------------------------------------------%
% edge-wise GWAS results
%--------------------------------------------------%
[pORA_eQTL_HCP, pORA_eQTL_Monash, pORA_eQTL_Consensus] = edgeGWAS_eQTL_ORA(); 