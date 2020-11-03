function plot_hubGroupsSurface(parc,nodeData,ts, side, hemi)
% plot hubs on surface
nodeLab = zeros(length(nodeData),1);

nodeLab(nodeData>ts(3)) = 1; 
nodeLab(nodeData>ts(2)) = 2; 
nodeLab(nodeData>ts(1)) = 3; 
modulesColor = GiveMeColors('degreeGroups'); 

boundColour = [82,82,82]/255; 

[f1, in] = plot_measureONsurface(parc, nodeLab, side , modulesColor, hemi);
hold on; 
%make boundaries on the surface
cdata_continuous = makeROIcont(in);
F = makeROIbound(cdata_continuous, in);
for i=1:length(F)
    plot3(F{i}(:,1), F{i}(:,2), F{i}(:,3), 'Color', boundColour, 'LineWidth',1.5)
end


end