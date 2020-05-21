function PlotVoronoiLandscapeC(C,P,clabel, markBestNets,yaxislabel)

mtype = {'sptl','neighbors','matching','clu-avg','clu-min','clu-max','clu-diff','clu-prod','deg-avg','deg-min','deg-max','deg-diff','deg-prod','com'};

if nargin < 4
   markBestNets = 0; 
end

if nargin < 5
    yaxislabel = 'Gamma';
end

    
    if iscell(C)
     for i = 1:13
         
         subplot(3,5,i)
         
        scatter(P{i}(:,1),P{i}(:,2),20,C{i},'filled')

if markBestNets
hold on
  [~,ORDIND] = sort(C{i},'descend');

scatter(P{i}(ORDIND(1:100),1),P{i}(ORDIND(1:100),2),5,'*k')

end
title(mtype{i})

ylabel(yaxislabel)
xlabel('Eta')
  c = colorbar;

c.Label.String = clabel;

colormap(turbo)
     end
    else
    
scatter(P(:,1),P(:,2),20,C,'filled')

if markBestNets
hold on
  [~,ORDIND] = sort(C,'descend');

scatter(P(ORDIND(1:100),1),P(ORDIND(1:100),2),5,'*k')

end
%title(model)

ylabel(yaxislabel)
xlabel('Eta')
  c = colorbar;

c.Label.String = clabel;

colormap(turbo)

    end

end



