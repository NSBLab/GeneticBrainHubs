% compare GCC for different link types
function geneScore = compareGCC(GCC, AdjMatrix,whatType, deg, hubThresh, doConnected)

if nargin<6
    doConnected=true; 
    fprintf('BY DEFAULT calculating scores on connected pairs of regions\n');
end
    
switch whatType
    case 'Connected'
        maskLinks1 = logical(AdjMatrix);
        maskLinks2 = logical(~AdjMatrix);
    case 'Unconnected'
        maskLinks1 = logical(~AdjMatrix);
        maskLinks2 = logical(AdjMatrix);
    case 'Rich' 
        isHub = deg>hubThresh; 
        maskHub = zeros(length(deg)); 
        maskHub(isHub, isHub) = 1;
        maskPeriph = zeros(length(deg)); 
        maskPeriph(~isHub, ~isHub) = 1;
        if doConnected
        maskLinks1 = logical(maskHub.*logical(AdjMatrix)); % connected rich links
        maskLinks2 = logical(maskPeriph.*logical(AdjMatrix)); % connected peripheral
        else
        maskLinks1 = logical(maskHub); % hub pairs
        maskLinks2 = logical(maskPeriph); % peripheral pairs
        end
    case 'RichFeeder'
        isHub = deg>hubThresh; 
        maskRF = zeros(length(deg)); 
        maskRF(isHub, isHub) = 1;
        maskRF(~isHub, isHub) = 1;
        maskRF(isHub, ~isHub) = 1;
        maskPeriph = zeros(length(deg)); 
        maskPeriph(~isHub, ~isHub) = 1;
        if doConnected
        maskLinks1 = logical(maskRF.*logical(AdjMatrix)); % connected rich and feeder links
        maskLinks2 = logical(maskPeriph.*logical(AdjMatrix)); % connected peripheral
        else
        maskLinks1 = logical(maskRF); % pairs of hubs and non-hubs
        maskLinks2 = logical(maskPeriph); % pairs of non-hubs
        end

        
end

geneScore = getGCCscore(GCC, maskLinks1, maskLinks2);
              
function scoreVal = getGCCscore(x, maskLinks1, maskLinks2)
    scoreVal = zeros(size(x,3),1); 
    for i=1:size(x,3)
        y = squeeze(x(:,:,i)); 
        y = masklHalf(y); 
        z1 = y(maskLinks1); % links that show an increase
        links1 = z1(~isnan(z1));
        z2 = y(maskLinks2); % links that do not show an increase
        links2 = z2(~isnan(z2));

        [~,~,~, stats] = ttest2(links1,links2, 'Vartype', 'unequal', 'Tail', 'right'); 
        scoreVal(i) = stats.tstat; 
    end
end
end