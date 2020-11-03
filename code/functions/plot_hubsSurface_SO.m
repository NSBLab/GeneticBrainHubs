function plot_hubsSurface_SO(parcellation,nodeData,ts, side, whatHemisphere)

nodeLab = zeros(length(nodeData),1);

nodeLab(nodeData>ts(3)) = 1; 
nodeLab(nodeData>ts(2)) = 2; 
nodeLab(nodeData>ts(1)) = 3; 
modulesColor = GiveMeColors('degreeGroups'); 

load('data/modules/FSAVERAGE_DATA_inflated.mat')
load('data/modules/fsaverage_surface_data.mat')

% Get the vertex ROI number

switch parcellation
    
    case 'random500'
        switch whatHemisphere
            case 'lh'
                parcdata = importVERTfile('data/modules/lh.random500.txt');
            case 'rh'
                parcdata = importVERTfile('data/modules/rh.random500.txt');
        end
        
    case 'random200'
        switch whatHemisphere
            case 'lh'
                parcdata = importVERTfile('data/modules/lh.random200.txt');
            case 'rh'
                parcdata = importVERTfile('data/modules/rh.random200.txt');   
        end

    case 'HCP'
        switch whatHemisphere
            case 'lh'
                parcdata = lh_HCPMMP1;
            case 'rh'
                parcdata = rh_HCPMMP1;    
        end

    case 'aparcaseg'
        switch whatHemisphere
            case 'lh'
                parcdata = lh_aparc;
            case 'rh'
                parcdata = rh_aparc;
        end
        
end

if strcmp(whatHemisphere, 'rh')
    surface.faces = rh_faces;
    surface.vertices = rh_inflated_verts; 
  
else
    surface.faces = lh_faces;
    surface.vertices = lh_inflated_verts; 
    
end

figure; 
set(gcf,'color','w');
set(gcf, 'Position', [500 500 950 750])

       
plotSurfaceROIBoundary(surface,parcdata,nodeLab,'midpoint',modulesColor,1,2);


if strcmp(side, 'inside') && strcmp(whatHemisphere, 'rh')
    turnAngle = -90;
elseif strcmp(side, 'inside') && strcmp(whatHemisphere, 'lh')
    turnAngle = 90;
elseif strcmp(side, 'outside') && strcmp(whatHemisphere, 'rh')
    turnAngle = 90;
elseif strcmp(side, 'outside') && strcmp(whatHemisphere, 'lh')
    turnAngle = -90;
end

view(turnAngle,0)
camlight('headlight')

axis off
axis tight
axis equal
end