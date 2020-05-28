function [f, out]  = plot_measureONsurface(parcellation, data, plotWhat, whatMap, whatHemisphere)
if nargin<4
    whatMap = 'magma';
    whatHemisphere = 'lh';
end

if nargin<5
    whatHemisphere = 'lh';
end

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
        
        
    case 'custom500'
        
     
        switch whatHemisphere
            case 'lh'
                parcdata = importVERTfile('data/modules/lh.custom500.txt');
            case 'rh'
                parcdata = importVERTfile('data/modules/rh.custom500.txt');
                
        end
        
    case 'cust250'
        
        %parcdata = custom500;
        switch whatHemisphere
            case 'lh'
                parcdata = importVERTfile('data/modules/lh.custom500.txt');
            case 'rh'
                parcdata = importVERTfile('data/modules/rh.custom500.txt');
                
        end
        
    case 'HCP'
        
        switch whatHemisphere
            case 'lh'
                parcdata = lh_HCPMMP1;
            case 'rh'
                parcdata = rh_HCPMMP1;
                
        end
        
        
    case 'custom200'
        
        switch whatHemisphere
            case 'lh'
                parcdata = lh_cust200;
            case 'rh'
                parcdata = rh_cust200;
                
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
    g.faces = rh_faces;
    g.vertices = rh_inflated_verts; 
  
else
    g.faces = lh_faces;
    g.vertices = lh_inflated_verts; 
    
end

gg = gifti;

% Map the data onto the ROI vertices
cdata = nan(size(parcdata));
for i = 1:length(data)
    cdata(parcdata == i) = data(i);
end

if strcmp(whatMap, 'magma')
    m = 200;
    cmap = [(magma(m)); 1 1 1];
else
    cmap = [(whatMap); 1 1 1];
    m = size(whatMap,1);
end


% For ROIs with no data give them a very small value so color mapping will work
%cdata(cdata == 0) = NaN;

cdata_nans = isnan(cdata);

cmin = nanmin(data);
cmax = nanmax(data);

cmapinds = linspace(cmin,cmax,m);
cmapspacing = cmapinds(m)-cmapinds(m-1);

cdata(cdata_nans) = cmax + cmapspacing;
gg.cdata = (cdata);

% Plot the zscored coefficients on the (inflated) fsaverage brain
% if strcmp(plotWhat, 'inside') && strcmp(whatHemisphere, 'rh')
% turnAngle = -90; 
% elseif strcmp(plotWhat, 'inside') && strcmp(whatHemisphere, 'lh')
% turnAngle = 90; 

if strcmp(plotWhat, 'inside') && strcmp(whatHemisphere, 'rh')
    turnAngle = -90;
elseif strcmp(plotWhat, 'inside') && strcmp(whatHemisphere, 'lh')
    turnAngle = 90;
elseif strcmp(plotWhat, 'outside') && strcmp(whatHemisphere, 'rh')
    turnAngle = 90;
elseif strcmp(plotWhat, 'outside') && strcmp(whatHemisphere, 'lh')
    turnAngle = -90;
end

        f=figure;
        set(gcf,'color','w');
        set(gcf, 'Position', [500 500 950 750])
        plot(g,gg)
        colormap(cmap)
        material dull
        view(turnAngle,0)
        camlight('headlight')

out.vertices = g.vertices;
out.faces = g.faces;

out.cdata = parcdata;

end