# GeneticBrainHubs
This repository provides Matlab code for reproducing figures from the paper entitled "Genetic influences on hub connectivity ofthe human connectome". The code was written using MATLAB_R2019b. Some functions are dependent on the version and should be accommodated by versions >R2017a.

Contact Aurina Arnatkeviciute by [email](mailto:aurina.arnatkeviciute@monash.edu).

### Data files
Data files required for this project are hosted on [this cloudstor repository](https://cloudstor.aarnet.edu.au/plus/s/yZafVFjELsQ8SoB).
Please download and unzip the data in the root directory.

### Reproducing figures
First, add all sub-folders to the path using `setupPaths()` function from the root directory. This will also create a `makeFigures` folder.
Code to reproduce the figures is located in the `generateFigures` folder. Script `generateFigures.m` calls all the functions to plot the figures. Run the scripts from the root directory. Figures will be saved in `makeFigures` folder in the root directory.


`Figure2()` - main text: heritability


`Figure3()` - main text: correlated gene expression


`Figure4()` - main text: generative modelling


`FigureS1()` - supplementary: topological properties of the connectome


`FigureS2_S5()` - supplementary: heritability-related analyses


`FigureS6_S11()` - supplementary: correlated gene expression related analyses


`FigureS12()` - supplementary: microstructure related analysis


`FigureS13()` - supplementary: generative modelling


`[pORA_eQTL_Monash, pORA_eQTL_HCP] = edgeGWAS_eQTL_ORA;` - edge-wise GWAS results
