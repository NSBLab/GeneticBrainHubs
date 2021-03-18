function plot_edges_brain(edge_values,nodeData, edge_cmap)

% load fsaverage annot file 
load('data/modules/FSAVERAGE_DATA_inflated.mat')
load('data/modules/fsaverage_surface_data.mat')

numNodes = size(edge_values,1); 
% 
surface_L.vertices = lh_pialverts2;
surface_L.faces = lh_faces;

surface_R.vertices = rh_pialverts2;
surface_R.faces = rh_faces;

NodesCoords = zeros(numNodes,3);

% Below is a simple boring way to get coordinates. Just take the average of
% the vertex coordinates which make up each ROI. Use if simpleton

% for i = 1:360/2
%    
%    NodesCoords(i,:) = mean(surface_L.vertices(lh_HCPMMP1==i,:));
%    NodesCoords(i+180,:) = mean(surface_R.vertices(rh_HCPMMP1==i,:));
%    
% end

% This is a fancier way. Within each ROI, find the vertex which is the
% closest to all other vertices. Use that vertex as the node coordinates

for i = 1:numNodes/2
   
   Roi_verts = surface_L.vertices(lh_HCPMMP1==i,:);
    
   [~,ROI_Vert_IND] = min(mean(squareform(pdist(Roi_verts))));
   
   NodesCoords(i,:) = Roi_verts(ROI_Vert_IND,:);
   
   Roi_verts = surface_R.vertices(rh_HCPMMP1==i,:);
    
   [~,ROI_Vert_IND] = min(mean(squareform(pdist(Roi_verts))));
   
   NodesCoords(i+numNodes/2,:) = Roi_verts(ROI_Vert_IND,:);
   
end
node_cmap = cbrewer('seq', 'Reds', 256); 
% select 500 edges at random
Plot3DConnectome(edge_values,NodesCoords,node_cmap,edge_cmap,nodeData,surface_L,surface_R);
% need to use camlight to make it pretty. This shines a light from the left
% and right

camlight('headlight')
% camlight(80,-10);
% camlight(-80,-10);

axis off
axis tight
axis equal
% Use view to change the position of the camera. view([90 0]) views the
% left hemisphere
% view([90 0])