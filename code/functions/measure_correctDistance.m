function [pairMeasCorr,FitCurve, c, data_uncorr] = measure_correctDistance(pairwiseMeasure, distanceMatrix, measureName, fitType, doPlot)
if nargin<4
    fitType = 'exp';
    doPlot = true; 
end

% this function performs distance correction on CGE data on a selected
% subset of genes, by regressing a bulk trend for all available genes

pairwiseMeasure = maskuHalf(pairwiseMeasure);
% get distances between regions
DISTmatrix = maskuHalf(distanceMatrix);

% get the fit
Fit{1} = fitType;
dataFIT = [DISTmatrix(:), pairwiseMeasure(:)];
[~,~,c] = GiveMeFit(DISTmatrix(:),pairwiseMeasure(:),Fit{1});
[FitCurveV] = getFitCurve(Fit,DISTmatrix(:),c);

% remove the fit from CGEmatrix on selected genes
Residuals = pairwiseMeasure(:) - FitCurveV;
FitCurve = reshape(FitCurveV,[length(pairwiseMeasure) length(pairwiseMeasure)]);
FitCurve(isnan(FitCurve)) = 0;
FitCurve = FitCurve+FitCurve';

% reshape into matrix for further use
pairMeasCorr = reshape(Residuals,[length(pairwiseMeasure) length(pairwiseMeasure)]);

% plot the distance relationship for visualisation
if doPlot
[~, y_uncor, ~, x_uncor] = BF_PlotQuantiles(DISTmatrix(:),pairwiseMeasure(:),26,1,1);   ...
     xlabel('Distance between regions (mm)'); ...
     ylabel(sprintf('%s', measureName));...
     set(gca,'fontsize',18);
     ylim([-1 1]);
hold on; scatter(DISTmatrix(:),FitCurveV,3, '.', 'r');

[~, y_cor, ~, x_cor]  = BF_PlotQuantiles(DISTmatrix(:),pairMeasCorr(:),26,1,1);   ...
     xlabel('Distance between regions (mm)'); ...
     ylabel(sprintf('%s', measureName));...
     set(gca,'fontsize',18);
     ylim([-1 1]);
     data_uncorr = table(x_uncor, y_uncor', x_cor, y_cor', 'VariableNames', {'Dist_uncorr','CGE_uncorr', 'Dist_corr','CGE_corr'}); 
end



end
