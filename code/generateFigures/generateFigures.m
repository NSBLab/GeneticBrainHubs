%--------------------------------------------------% 
% This script runs all the functions to reproduce figures in the manuscript
%--------------------------------------------------%

setupPaths()
%--------------------------------------------------%
% Figure 2 - main text: heritability and supplementary figure S2 visualising spatial distribution of most and least heritable links. 
%--------------------------------------------------%
Figure2_and_FigureS2()

%--------------------------------------------------%
% Figure 3 - main text: correlated gene expression
%--------------------------------------------------%
Figure3()

%--------------------------------------------------%
% Figure 4 - main text: generative modelling (space and topology) when optimising max KS.
%--------------------------------------------------%
Figure4('energy')

%--------------------------------------------------%
% Figure 5 - main text: generative modelling (space, topology and genes). 
%--------------------------------------------------%
Figure5()

%--------------------------------------------------%
% Figure S1 - supplementary: topological properties 
% of the connectome
%--------------------------------------------------%
FigureS1()

%--------------------------------------------------%
% Figure S3-S5 - supplementary: heritability-related analyses
%--------------------------------------------------%
FigureS3_S5()

%--------------------------------------------------%
% Figure S6-S11 - supplementary: correlated gene 
% expression related analyses
%--------------------------------------------------%
FigureS6_uncorrectedCGE()
FigureS6_S11()

%--------------------------------------------------%
% Figure S12 - supplementary: generative modelling (space and topology) 
% when optimising degree correlation.
%--------------------------------------------------%
FigureS12()

%--------------------------------------------------%
% Figure S13 - supplementary: heritability for rich, feeder and peripheal links 
% where the significance is evaluated using label randomisation.
%--------------------------------------------------%
FigureS13()

%--------------------------------------------------%
% Figure S14 - supplementary: genetic variance.
%--------------------------------------------------%
FigureS14()
