function getMaxVal = RichClubMouse(Adj,pairwiseMeasure,nodeData, whatTail, whatDistribution, colorOut, colorIn)
% ------------------------------------------------------------------------------
% Function plots coexpression for rich/feeder/peripheral lins as a function
% of degree using mean to summarise coexpression at each threshold
%-------------------------------------------------------------------------------
% INPUTS:
% ------------------------------------------------------------------------------

realLinkData = pairwiseMeasure;
numBins = 'all'; % Range of k to plot the rich club curve across

% ------------------------------------------------------------------------------
% Assign data measured at each link in the network
% ------------------------------------------------------------------------------
linkedAdj = Adj;
% Add a mask to only include particular types of connections
numNeurons = size(linkedAdj,1);
nn = linspace(1,numNeurons,numNeurons);
allLinkData = realLinkData;
allLinkData(Adj==0) = NaN;
%allLinkData(allLinkData<0.0001) = NaN; 

numNodes = length(Adj);

% ------------------------------------------------------------------------------
% Get groups of links based on their degree (or use all k)
% ------------------------------------------------------------------------------
sortK = sort(nodeData,'descend');
maxK = sortK(2); % Up to the second-highest k
if strcmp(numBins,'all')
    % All k are a bin
    kr = min(nodeData):maxK;
else % numeric number of bins
    kr = linspace(min(nodeData),maxK,numBins);
    kr = round(kr);
end
krAll = min(nodeData):max(nodeData);

% ------------------------------------------------------------------------------
% Go through each class of links and compute statistics on the set of link
% data compared to the nulls:

whatLinks = {'rich','feeder','local'};
% Ben Fulcher, 2015-01-06
% Computes a t-test between special and non-special links for each k, and each link-type:
tStats = zeros(length(kr),length(whatLinks));
allHubHub = cell(1,1); % make a 1-component cell for consistency with null version
allHubHub{1} = cell(length(kr),length(whatLinks));
for i = 1:length(kr)
    for j = 1:length(whatLinks) % loop across rich, feedin, feedout, and local connections:

        % Keep links between high degree nodes that have linkData (i.e., links that exist):
        % (assumption that nodeData will be the same for all null networks)
        % Ben Fulcher, 2015-01-06 -- changed from >= to > to match actual rich-club coeff definition
        r = (nodeData > kr(i));

        keepMe = logical(zeros(numNodes,numNodes));
        % Add ones for a given type of link:
        switch j
            case 1 % 'rich'
                keepMe(r,r) = 1;
            case 2 % 'feeder'
                keepMe(~r,r) = 1;
                keepMe(r,~r) = 1;
            case 3 % 'local'
                keepMe(~r,~r) = 1;

        end
        % Remove missing data (NaNs encode no link):
        keepMe(isnan(allLinkData)) = 0;

        % Keep values assigned to each remaining link as this element of hubhubData
        linkDataSpecial = allLinkData(keepMe);
        allHubHub{1}{i,j} = linkDataSpecial;

        notSpecial = (~isnan(allLinkData) & ~keepMe);
        linkDataNotSpecial = allLinkData(notSpecial);

        % 2-sample t-test for special links greater than non-special links:

        [~,p] = ttest2(linkDataSpecial,linkDataNotSpecial,'Vartype','unequal', 'Tail',whatTail);

        tStats(i,j) = p;
    end
end

% ------------------------------------------------------------------------------
% Plot as rich plots
% ------------------------------------------------------------------------------
myColors = GiveMeColors('RFPU'); 
% % reorder colours for correct overlay, so rich are plotted at the top
% % [BF_getcmap('spectral',4,1),BF_getcmap('set2',4,1)];
% myColors(1,:) = myColorsOrig(3,:); 
% myColors(2,:) = myColorsOrig(2,:); 
% myColors(3,:) = myColorsOrig(1,:); 

plotOnOne = true; % plot all on one figure
includeHist = true;
plotJustRich = false;
sigThresh = 0.05;

for j = 1:length(whatLinks)
    if plotJustRich && (j < 1)
        break
    end

        if includeHist
            if (~plotOnOne || (plotOnOne==1 && j==1))


                figure('color','w');
                sp=subplot(5,3,1:6);

                switch whatDistribution
                    case 'barCount'

                        N = arrayfun(@(x)sum(nodeData==x),krAll);
                        bar(krAll,N,'EdgeColor',colorOut,'FaceColor',colorOut);

                    case 'histogram'

                         histogram(nodeData,35,'EdgeColor',[1 1 1],'FaceColor',colorOut, 'FaceAlpha',1);

                    case 'kernel'

                        pdSix = fitdist(nodeData','Kernel','BandWidth',3);
                        x = 0:1:max(nodeData)+5;
                        ySix = pdf(pdSix,x);
                        plot_shaded(x,ySix, 'Color', colorOut)
                        %plot(x,ySix,'k-','LineWidth',2)
                end

                xlim([min(nodeData)-0.8,max(nodeData)+0.8]);
                xticks([]); box off;
                ylabel('Frequency', 'FontSize', 18)
                get(gca, 'YTick');
                set(gca, 'FontSize', 18)
                sp=subplot(5,3,7:15);
                hold on;

            end
        else
            if (~plotOnOne || (plotOnOne==1 && j==1))
                figure('color','w'); hold on
            end
        end



    % Plot flat line for rich
    % rich are now labeled 3
    if j==1

        plot([kr(1),krAll(end)],ones(2,1)*nanmean(allHubHub{1}{1,j}),':','color','k','LineWidth',3)

    end

    valsPlot = allHubHub{1}(:,j);
    realTrajectory = cellfun(@nanmean,valsPlot);
    xlim([min(nodeData)-0.5,max(nodeData)+0.5]);


    % p-values from 2-sample t-test with unequal variances:
    pvalues = tStats(:,j);

    isSig = (pvalues < sigThresh); % significantly higher than null

    % meadian (real data trajectory):
    markerStyle = 'o';

    % meadian trajectory:
    % +/- std:
    % find max number of values across all distributions for a selected
    % link type
    Lv = zeros(length(valsPlot),1); 
    for w=1:length(valsPlot)
        Lv(w) = length(valsPlot{w});
    end

    A = nan(length(kr), max(Lv));
    for p=1:length(kr)
        V = valsPlot{p};
        for q=1:length(V)

            A(p,q) = valsPlot{p}(q);

        end
    end
    
    % NaNs can't be plotted, in that case make x and y variables shorter
    INDkeep = find(~isnan(nanmean(A,2))); 
    if ~isempty(INDkeep)
    plot_distribution(kr(INDkeep),A(INDkeep,:)', 'Color', myColors(j,:)); 
    hold on;
    end

    if any(isSig)
        plot(kr(isSig),realTrajectory(isSig),markerStyle,'MarkerEdgeColor',myColors(j,:),...
            'MarkerFaceColor',colorIn,'LineWidth',1.5,'MarkerSize',6)
    end

    xLimits = get(gca,'xlim'); yLimits = get(gca,'ylim');

%     if ~plotJustRich
%         text(xLimits(1)+0.1*diff(xLimits),yLimits(1)+0.9*diff(yLimits)-j/20,whatLinks{j},'color',myColors(j,:),'FontSize',18)
%     end


end

set(gcf, 'Position', [500 500 750 500])

axisName = {'Mean correlated', 'gene expression'};
ylabel(axisName, 'FontSize', 18)
xlabel('Node degree, k','FontSize', 18);
getMaxVal = cellfun(@nanmean,allHubHub{1}(:,1));
getMinVal = cellfun(@nanmean,allHubHub{1}(:,3));
ylim([min(getMinVal) max(getMaxVal)+nanstd(getMaxVal)])
set(gca,'fontsize', 18);
end
