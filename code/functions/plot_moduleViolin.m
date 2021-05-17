function [res,f, S_export] = plot_moduleViolin(allData, measureModRICH, netNamesAll, netassignments, pairwiseMeasure, nodeData, groupAdjlog, plotWhat, khub)

% colours for modules
c = GiveMeColors('RFPU');
c = c(1:3,:);
c = [c;[1 1 1]; [1 1 1]; c];

% [1 1 1] colours are for gaps

modulesColor = GiveMeColors('funcModules'); 
modulesColorALL = vertcat(modulesColor, [1 1 1], [1 1 1], c);

switch plotWhat
    case 'CGE'
        adj = groupAdjlog(1:180, 1:180);
        nodeData = nodeData(1:180);
        
        matrix = pairwiseMeasure;
        yName = {'Spatially-corrected', 'correlated gene expression'};
        yRange = [-1 1];
    case 'heritability'
        adj = groupAdjlog;
        % this will be determined outside of the script
        %pairwiseMeasure(pairwiseMeasure==0) = NaN;
        matrix = pairwiseMeasure;
        yName = 'Heritability';
        yRange = [0 1];
end

% make violin plots for inter-moduler and intra-modular rich links
[MeasureINT, res, MeasureINTnames] = compareRich_InterIntra_modules(matrix, netassignments, adj, nodeData, khub);


% here combine modules and inter/intra modular data into one cell
allData = vertcat(allData{1:length(modulesColor)}, 2, 2, MeasureINT(1:3)', 2, 2, MeasureINT(4:6)');
netNamesAll = vertcat(netNamesAll{1:length(modulesColor)},{'gap11'}, {'gap12'}, MeasureINTnames(1:3)',{'gap21'}, {'gap22'}, MeasureINTnames(4:6)'); 

f=figure('color','white');
set(gcf, 'Position', [10 10 1200 500]);

for i=1:length(allData)
    switch plotWhat
        case 'heritability'
            moduleVals = allData{i};
            % this will be determined outside of the function as an input
            %moduleVals(moduleVals==0) = NaN;
            S.(netNamesAll{i}) = moduleVals;
        case 'CGE'
            S.(netNamesAll{i}) = allData{i};
    end
 
end

for i=1:length(modulesColor)

    if ~isempty(measureModRICH{i})
        Sr.(netNamesAll{i}) = measureModRICH{i};
    else
        Sr.(netNamesAll{i}) = 0;
    end
    
end

for i=length(modulesColor):length(allData)
    
    Sr.(netNamesAll{i}) = 0;
    
end

set(gca,'FontSize',16)
violins = violinplot(S);
set(gcf,'color','w');

ylabel(yName)
ylim(yRange)
% change colors

for i=1:size(modulesColorALL,1)
    violins(1,i).ViolinColor = modulesColorALL(i,:);
    violins(1,i).EdgeColor = modulesColorALL(i,:);
    violins(1,i).BoxColor = [.25 .25 .25];
end

hold on;
% add rich links in infferent colour
violinsR = violinplot(Sr);
% change colors to dark grey
for i=1:size(modulesColorALL,1)
    % find x coordinates for dots
    %if isempty(violinsR(1,i).ScatterPlot) && isempty(violinsR(1,i).MeanPlot)
    
    if isempty(violinsR(1,i).ScatterPlot) && ~violinsR(1,i).MedianPlot.YData==0

        violinsR(1,i).MedianPlot.Visible = 'on';
        
        
        violinsR(1,i).MedianPlot.MarkerFaceColor = [0.25 0.25 0.25];
        violinsR(1,i).MedianPlot.MarkerEdgeColor = [0.25 0.25 0.25];
        violinsR(1,i).MedianPlot.SizeData = 23; 
        violinsR(1,i).MedianPlot.MarkerFaceAlpha = 1;
        
    elseif isempty(violinsR(1,i).ScatterPlot) && violinsR(1,i).MedianPlot.YData==0
        
        violinsR(1,i).MedianPlot.Visible = 'off';
        
    elseif ~isempty(violinsR(1,i).ScatterPlot)
        
        [~,indALL] = ismember(violinsR(1,i).ScatterPlot.YData,violins(1,i).ScatterPlot.YData);
        %[~, indALL] = intersect(violins(1,i).ScatterPlot.YData, violinsR(1,i).ScatterPlot.YData, 'stable');
        % get corresponding X values from the main plot
        Xvals = violins(1,i).ScatterPlot.XData(indALL);
        % use those values for a rich plot
        violinsR(1,i).ScatterPlot.XData = Xvals;
        violinsR(1,i).MedianPlot.Visible = 'off';
        violinsR(1,i).ScatterPlot.MarkerFaceAlpha = 1;
        % change dot colour
        violinsR(1,i).ScatterPlot.MarkerFaceColor = [0.25 0.25 0.25];
        
    end
    % remove violin
    violinsR(1,i).ViolinPlot.Faces=NaN;
    violinsR(1,i).WhiskerPlot.Visible='off';
    violinsR(1,i).BoxPlot.Visible = 'off';
    
end

% add number of points to tick label

%INTlabels ={'intra-rich','intra-feeder', 'intra-peripheral', 'gap21', 'gap22', 'inter-rich', 'inter-feeder','inter-peripheral'}; 
INTlabels ={'rich','feeder', 'peripheral', 'gap21', 'gap22', 'rich', 'feeder','peripheral'}; 
netNamesAll = vertcat(netNamesAll{1:length(modulesColor)+2},INTlabels'); 

for l=1:size(allData,1)
    
    TF = startsWith(netNamesAll{l},'gap'); 
    if TF
        A = {''}; 
    else  
        A = strcat(netNamesAll{l}, {' '}, {'('}, num2str(length(allData{l})),{')'}); 
    end
    violinLabels{l} = A{1};

end

ax = gca; 
% remobe ticks from only X axis
ax.XAxis.TickLength = [0,0];
xticklabels(violinLabels); 
xtickangle(30)

S_export = rmfield(S,{'gap11', 'gap12', 'gap21', 'gap22'});

end
