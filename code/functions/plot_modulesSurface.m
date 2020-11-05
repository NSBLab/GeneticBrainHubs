function [f1, f2, f3, f4] = plot_modulesSurface(netassignments)
% plot modules on surface

modulesColor = GiveMeColors('funcModules'); 

f1 = plot_hubsSurface_modules_SO('HCP', netassignments(1:180), modulesColor, 'outside', 'lh');

f2 = plot_hubsSurface_modules_SO('HCP', netassignments(181:360), modulesColor,'inside', 'rh');

f3 = plot_hubsSurface_modules_SO('HCP', netassignments(1:180), modulesColor, 'inside', 'lh');

f4 = plot_hubsSurface_modules_SO('HCP', netassignments(181:360), modulesColor, 'outside', 'rh');

end