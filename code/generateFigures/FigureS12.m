%--------------------------------------------------% 
% Figure S12
%--------------------------------------------------%
function FigureS12()

weight = 'standard'; %for GenCog'standard'; 
parc = 'random500';
op = selectCONmetrics(parc, weight); 


plotOptions.colIn = [1 1 1];
plotOptions.colorOut = [82 82 82]/255; %[0,69,41]/255; %[.35 .35 .35]; %[4 90 141]/255; %[.45 .45 .45]; %[5 113 176]/255; 
plotOptions.whatDistribution = 'histogram'; 

[coexpData, A, matrices, coordinates, avWeight] = giveConnExp_HCP(parc,op.tract,op.probe,weight,op.brainPart,op.nRem);
[~, ~, M, ~, avWeightFA] = giveConnExp_HCP(parc,op.tract,op.probe,'FA',op.brainPart, 0);
    
numLC = size(coexpData.averageCoexpression,1); 

[GrSC, mDIST] = giveMeRichClub(matrices, coordinates, op.groupConn ,op.densThreshold, false, op.cvMeasure, op.consThr);
groupAdjlog = logical(GrSC); 
nodeData = degrees_und(groupAdjlog);

% Load BigBrain data
yregress = true; 
[BBmean,microDist,distNorm, BBmpc, mp, BBskew] = getBBdata('random500', yregress);

RichClubHuman(groupAdjlog,BBmpc,nodeData, 'right', plotOptions.whatDistribution, plotOptions.colorOut, plotOptions.colIn);
ylabel('Mean microstructural covariance')
set(gcf, 'Position', [500 500 750 550])
set(gca,'fontsize', 18);
ylim([0.6 1.1])

figureName = sprintf('makeFigures/MPCmean_%s_%d.png', parc, round(op.densThreshold*100));
print(gcf,figureName,'-dpng','-r600');
end
