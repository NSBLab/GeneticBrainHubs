function plot_expressionSurface(parc,nodeData,side, hemi)
% plot gene expression on surface
modulesColor = cbrewer('div', 'RdBu', 25);

boundColour = [82,82,82]/255; 

[f1, in] = plot_measureONsurface(parc, nodeData, side, modulesColor, hemi);
hold on; 
%make boundaries on the surface
cdata_continuous = makeROIcont(in);
F = makeROIbound(cdata_continuous, in);
for i=1:length(F)
    plot3(F{i}(:,1), F{i}(:,2), F{i}(:,3), 'Color', boundColour, 'LineWidth',0.5)
end


end