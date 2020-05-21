function [f1, f2, f3, f4] = plot_modulesSurface(netassignments)
% plot modules on surface

modulesColor = GiveMeColors('funcModules'); 

[f1, in] = plot_measureONsurface('HCP', netassignments(1:180), 'outside', modulesColor, 'lh');
hold on;
% make boundaries on the surface
cdata_continuous = makeROIcont(in);
F = makeROIbound(cdata_continuous, in);
for i=1:length(F)
    plot3(F{i}(:,1), F{i}(:,2), F{i}(:,3), 'k','LineWidth',2)
end


[f2, in] = plot_measureONsurface('HCP', netassignments(181:360), 'inside', modulesColor, 'rh');
cdata_continuous = makeROIcont(in);
F = makeROIbound(cdata_continuous, in);
hold on;

for i=1:length(F)
    plot3(F{i}(:,1), F{i}(:,2), F{i}(:,3), 'k','LineWidth',2)
end

[f3, in] = plot_measureONsurface('HCP', netassignments(1:180), 'inside', modulesColor, 'lh');
cdata_continuous = makeROIcont(in);
F = makeROIbound(cdata_continuous, in);
hold on;

for i=1:length(F)
    plot3(F{i}(:,1), F{i}(:,2), F{i}(:,3), 'k','LineWidth',2)
end

[f4, in] = plot_measureONsurface('HCP', netassignments(181:360), 'outside', modulesColor, 'rh');
cdata_continuous = makeROIcont(in);
F = makeROIbound(cdata_continuous, in);
hold on;

for i=1:length(F)
    plot3(F{i}(:,1), F{i}(:,2), F{i}(:,3), 'k','LineWidth',2)
end
end