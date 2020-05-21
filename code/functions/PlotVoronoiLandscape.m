function PlotVoronoiLandscape(E,P,markBestNets,yaxislabel)


mtype = {'sptl','neighbors','matching','clu-avg','clu-min','clu-max','clu-diff','clu-prod','deg-avg','deg-min','deg-max','deg-diff','deg-prod','com'};

if nargin < 3
   markBestNets = 0; 
end

if nargin < 4
    yaxislabel = 'Gamma';
end

    crange = [0 1];


    clabel = 'Energy';

    if iscell(E)
        
        
        for i = 1:13
            
            subplot(3,5,i)
    
scatter(P{i}(:,1),P{i}(:,2),20,E{i},'filled')

if markBestNets == 1
hold on
  [~,ORDIND] = sort(E{i},'ascend');

scatter(P{i}(ORDIND(1:100),1),P{i}(ORDIND(1:100),2),5,'*k')

end

title(mtype{i})

ylabel(yaxislabel)
xlabel('Eta')
  c = colorbar;

c.Label.String = clabel;
if ~isempty(crange)
caxis(crange)
end
colormap(turbo)

        end

    else
    
scatter(P(:,1),P(:,2),20,E,'filled')

if markBestNets == 1
hold on
  [~,ORDIND] = sort(E,'ascend');

scatter(P(ORDIND(1:100),1),P(ORDIND(1:100),2),5,'*k')

end



ylabel(yaxislabel)
xlabel('Eta')
  c = colorbar;

c.Label.String = clabel;
if ~isempty(crange)
caxis(crange)
end
colormap(turbo)

    end
end



