%--------------------------------------------------%
% Figure S1
%--------------------------------------------------%
function FigureS1()

whatDWI = 'HCP';
weight = 'standard';
parcellation = 'HCP';

op = selectCONmetrics(parcellation, weight);

if strcmp(whatDWI, 'HCP')
    
    [coexpData, A, matrices, coordinates, avWeight] = giveConnExp_HCP(parcellation,op.tract,op.probe,weight,op.brainPart,op.nRem);
    [~, ~, M, ~, avWeightFA] = giveConnExp_HCP(parcellation,op.tract,op.probe,'FA',op.brainPart, 0);
    
elseif strcmp(whatDWI, 'GenCog')
    
    [coexpData, A, matrices, coordinates, avWeight] = giveConnExp_GenCog(parcellation,op.tract,op.probe,weight,op.brainPart,op.nRem);
    [~, ~, M, ~, avWeightFA] = giveConnExp_GenCog(parcellation,op.tract,op.probe,'FA',op.brainPart, 0);
    
end

yVals = [0.95 1.35];
yValsFA = [0.95 1.15];
numRepeats = 1000;
numShuffle = 50;
whatDistribution = 'histogram';

% choose colors
colorOut = [82 82 82]/255; % outside of the circles
colorIn = [1 1 1]; % inside of the circles

[GrSC, mDIST] = giveMeRichClub(matrices, coordinates, op.groupConn, op.densThreshold, false, op.cvMeasure, op.consThr, 50, 100 , 'bu', 'randmio_und', yVals);
GrFA = avWeightFA.*logical(GrSC);
GrSClog = log(GrSC);
GrSClog(isinf(GrSClog)) = 0;
nodeDeg = degrees_und(GrSC);

% a) topological RC
[~, dMiddle_top, PhiNormMean_top] = PlotRichClub(GrSC,mDIST,'bu','randmio_und', numShuffle, numRepeats, yVals, whatDistribution, colorOut, colorIn);
figureName = sprintf('makeFigures/RCbin_%s_%d.png', parcellation, round(op.densThreshold*100));
print(gcf,figureName,'-dpng','-r300');

% b) weighted RC - log(SC) weight, topology fixed, weights randomised
[~, dMiddle_wsc, PhiNormMean_wsc] = PlotRichClub(GrSClog,mDIST,'wu','shuffleWeights', numShuffle, numRepeats, yVals, whatDistribution, colorOut, colorIn);
figureName = sprintf('makeFigures/RCwei_logSC_%s_%d.png', parcellation, round(op.densThreshold*100));
print(gcf,figureName,'-dpng','-r300');

% c) weighted RC - FA weight, topology fixed, weights randomised
[~, dMiddle_wfa, PhiNormMean_wfa] = PlotRichClub(GrFA,mDIST,'wu','shuffleWeights', numShuffle, numRepeats, yValsFA, whatDistribution, colorOut, colorIn);
figureName = sprintf('makeFigures/RCwei_FA_%s_%d.png', parcellation, round(op.densThreshold*100));
print(gcf,figureName,'-dpng','-r300');

% d) connection length as a function of degree
RichClubHuman_TOPO(GrFA,mDIST,nodeDeg, true, whatDistribution, colorOut, colorIn);
axisName = {'Mean connection', 'distance (mm)'};
ylabel(axisName, 'FontSize', 18)
xlabel('Node degree, k','FontSize', 18);
ylim([50 65])
figureName = sprintf('makeFigures/DIST_%s_%d.png', parcellation, round(op.densThreshold*100));
print(gcf,figureName,'-dpng','-r300');

% e) binary normalised edge communicability as a function of degree
% - null, is randomising topology
ComNormBIN = normCommunicability(GrFA, 'bin', numShuffle, numRepeats, 'randmio_und');
RichClubHuman_TOPO(GrFA,ComNormBIN,nodeDeg, true, whatDistribution, colorOut, colorIn);
axisName = {'Mean normalised binary', 'connection communicability'};
ylabel(axisName, 'FontSize', 18)
xlabel('Node degree, k','FontSize', 18);
figureName = sprintf('makeFigures/COMMbin_%s_%d.png', parcellation, round(op.densThreshold*100));
print(gcf,figureName,'-dpng','-r300');

% f) weighted normalised edge communicability as a function of degree
% - null, is randomising weights while keeping the topology
ComNormWEI = normCommunicability(GrFA, 'wei', numShuffle, numRepeats, 'shuffleWeights');
RichClubHuman_TOPO(GrFA,ComNormWEI, nodeDeg, true, whatDistribution, colorOut, colorIn);
axisName = {'Mean normalised weighted', 'connection communicability'};
ylabel(axisName, 'FontSize', 18)
xlabel('Node degree, k','FontSize', 18);
figureName = sprintf('makeFigures/COMMwei_%s_%d.png', parcellation, round(op.densThreshold*100));
print(gcf,figureName,'-dpng','-r300');
end



