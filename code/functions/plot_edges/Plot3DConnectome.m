function p = Plot3DConnectome(A,node_coords,NodeCmap,EdgeCmap,Nodedata,surface_L,surface_R,NodeCmapLimits,EdgeCmapLimits)

% Input:
%       A = N*N Adjacency matrix
%
%       node_coords = N*3 matrix of node coordinates
%
%       NodeCmap = a colourmap to apply to nodes or and RGB values to use
%       for all nodes
%
%       EdgeCmap = a colourmap to apply to edges or and RGB values to use
%       for all edges. The colour of each edge is determined by its value
%       in A
%
%       NodeCdata = A vector of length N specifying the colour data of 
%       nodes (will map onto NodeCmap)
%
%       surface_L = a structure containing the fields 'vertices' and
%       'faces'. 'vertices' is the vertices of the left surface mesh and
%       'faces' are the faces of the left surface mesh
%
%       surface_R = a structure containing the fields 'vertices' and
%       'faces'. 'vertices' is the vertices of the right surface mesh and
%       'faces' are the faces of the right surface mesh
%
%       NodeCmapLimits = user specified limits of NodeCmap (i.e., the range
%       of values you want the colourmap to apply to) 
%
%       EdgeCmapLimits = user specified limits of EdgeCmap (i.e., the range
%       of values you want the colourmap to apply to) 

% Stuart Oldham, Monash University, 2021

% Number of nodes
NNodes = length(A);

if ~exist('Nodedata','var') || isempty(Nodedata)
    
    Nodedata = ones(NNodes,1)*0.75;
    
end

% Set the limits of the node colormap

if ~exist('NodeCmapLimits','var')

    NodeCmapLimits = [min(Nodedata) max(Nodedata)];

end

% Find the RGB value of each node

NodeRGB = vals2colormap(Nodedata, NodeCmap, NodeCmapLimits);

% Set up Node data

[x,y,z] = sphere(64);
fvc = surf2patch(x,y,z,'triangles');

% Define the size of all nodes. Can easily change the code to have a 

NodeSize = Nodedata*5/max(Nodedata); 

% Define the radius of each edge

EdgeRadius = 0.15;

%% Start the plotting!

% Plot the left hemisphere if it is available

if exist('surface_L','var')

    p = patch(surface_L);
    set(p,'EdgeColor','none','FaceColor',[.9 .9 .9],'FaceAlpha',.25);
    material dull
    hold on

end

% Plot the right hemisphere if it is available

if exist('surface_R','var')
    
    p2 = patch(surface_R);
    set(p2,'EdgeColor','none','FaceColor',[.9 .9 .9],'FaceAlpha',.25);
    material dull
    hold on

end

% Plot the nodes

NODES.axes1 = gca;
for j = 1:NNodes
    fvc_temp = fvc;
    fvc_temp.vertices = fvc_temp.vertices.*NodeSize(j) + node_coords(j,:); 
    NODES.s(j) = patch(fvc_temp,'EdgeColor','none','FaceColor',NodeRGB(j,:),'Clipping','off');
    material dull
    hold on
end

% Get the upper triangle values

[A_vec,UPPER_INDS] = triu2vec(A,1);

[x,y] = ind2sub([NNodes NNodes],UPPER_INDS(A_vec>0));

EdgeVals = A_vec(A_vec>0);

% Set the limits of the edge colormap. Colormap is only applied to edges
% that exist

if ~exist('EdgeCmapLimits','var')

    EdgeCmapLimits = [min(EdgeVals) max(EdgeVals)];

end

EdgeRGB = vals2colormap(EdgeVals, EdgeCmap, EdgeCmapLimits);

% Plot the edges

EDGES.axes1 = gca;
for i = 1:length(x)

    NODE1_COORDS = node_coords(x(i),:);
    NODE2_COORDS = node_coords(y(i),:);
    
    % Make a cylinder where one side starts at one node location, and the
    % other side ends at th eother node location
    
    [X, Y, Z] = cylinder2P(EdgeRadius, 100,NODE1_COORDS,NODE2_COORDS);
    fvc_edges = surf2patch(X,Y,Z,'triangles');
    
    EDGES.edge(i) = patch(fvc_edges,'FaceColor',EdgeRGB(i,:),'EdgeColor','none','Clipping','off');

end

material dull

