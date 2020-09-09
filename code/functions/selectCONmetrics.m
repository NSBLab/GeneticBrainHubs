
function op = selectCONmetrics(parc, weight)
% selectCONmetrics    Defines parameters for connectivity and gene 
% expression data
%
%---INPUTS:
% parc, the name of the cortical parcellation used
% weight, type of connection weights in the connectome - FA or SC
%
%---OUTPUT:
% op, structure containing default parameters

% ------------------------------------------------------------------------------

% ------------------------------------------------------------------------------
switch parc
        
    case 'HCP'
        densThreshold = 0.2;
        khub = 105;
    case 'Schaefer400'
        densThreshold = 0.2;
        khub = 105;
        
    case 'random500'
        densThreshold = 0.1;
        khub = 70;
end

switch weight
    case 'SC'
        weight = 'standard';
    case 'FA'
        weight = 'FA';
end

op.densThreshold = densThreshold; 
op.groupConn = 'CVmeasure'; 
op.cvMeasure = 'strength'; 
op.calcRC = false; 
op.consThr = 0.3; 
op.tract = 'iFOD2'; 
op.probe = 'RNAseq'; 
op.conW = 'standard'; 
op.brainPart = 'LRcortex'; 
op.nRem = 0; 
op.weight = weight; 
op.khub = khub; 

end