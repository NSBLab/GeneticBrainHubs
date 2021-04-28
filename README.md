# GeneticBrainHubs

![image](https://user-images.githubusercontent.com/14971156/116356515-b8fd9f80-a83e-11eb-8e7a-5dcfe3f7ba79.png)

This repository provides Matlab code for reproducing figures from the paper entitled:

- Arnatkeviciute et al. (2021) [:green_book: 'Genetic influences on hub connectivity ofthe human connectome'](https://doi.org/10.1101/2020.06.21.163915).

The code was written using MATLAB_R2019b. Some functions are dependent on the version and should be accommodated by versions >R2017a. For earlier versions, please replace all `"` with `'` in `importHeritabilityResultwP.m` and `importGENEIDfile.m` functions. 

Dependencies: [Brain Connectivity Toolbox](https://sites.google.com/site/bctnet/), Version 2017-15-01. Included in this repository.

Contact Aurina Arnatkeviciute by [email](mailto:aurina.arnatkeviciute@monash.edu).

### Data files
Data files required for this project are hosted in on the associated [Zenodo repository](https://doi.org/10.5281/zenodo.4724186).
Please download and unzip the data in the root directory.

### Reproducing figures
First, add all sub-folders to the path using `setupPaths()` function from the root directory. This will also create a `makeFigures` folder.
Code to reproduce the figures is located in the `generateFigures` folder. Script `generateFigures.m` calls all the functions to plot the figures. Run the scripts from the root directory. Figures will be saved in `makeFigures` folder in the root directory.


`Figure2_and_FigureS2()` - main text: heritability and supplementary figure S2 visualising spatial distribution of most and least heritable links. 

`Figure3()` - main text: correlated gene expression.

`Figure4('energy')` - main text: generative modelling (space and topology) when optimising max KS.

`Figure5()` - main text: generative modelling (space, topology and genes). 

`FigureS1()` - supplementary: topological properties of the connectome.

`FigureS3_S5()` - supplementary: heritability-related analyses. 

`FigureS6_uncorrectedCGE()` - supplementary: correlated gene expression without correction for distance effects.

`FigureS6_S11()` - supplementary: correlated gene expression related analyses. 

`FigureS12()` - supplementary: generative modelling (space and topology) when optimising degree correlation.

`FigureS13()` - supplementary: heritability for rich, feeder and peripheal links where the significance is evaluated using label randomisation. 

`FigureS14()` - supplementary: genetic variance.


