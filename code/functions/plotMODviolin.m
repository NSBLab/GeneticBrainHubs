% make violin plot for modelling results
% make a violin plot
function [S,I] = plotMODviolin(ENERGY, num, mtype, selWhat, doPlot)
if nargin<5
    doPlot = true; 
end

S = struct;  
% assign colours based on original order
if length(ENERGY) == 13
    ModelColors = zeros(length(mtype),3);
    ModelColors(1,:) = [44,123,182]/255;
    for i = 2:3
        ModelColors(i,:) = [215,25,28]/255;
    end
    for i = 4:8
        ModelColors(i,:) = [253,174,97]/255;
    end
    for i = 9:13
        ModelColors(i,:) = [103,169,207]/255;
    end
    ModelColors(14,:) = [255 192 0]/255;
    
    % remove last, there are only 13 models
    ModelColors(14,:) = [];
elseif length(ENERGY) == 5
    
    ModelColors = zeros(length(mtype),3);
    ModelColors(1,:) = [103,169,207]/255; % ST, same as degree average colour
    ModelColors(2,:) = [90,174,97]/255; % models with genes - green
    ModelColors(3,:) = [90,174,97]/255; % models with genes
    ModelColors(4,:) = [44,123,182]/255; % spatial
    ModelColors(5,:) = [90,174,97]/255; % models with genes

end


% order models based on mean fit values
for j=1:size(ENERGY,2)
    switch selWhat
        case 'highest'
            [~,indM] = sort(ENERGY{j},'descend');
        otherwise    
            [~,indM] = sort(ENERGY{j});
    end
     
    lowIND1 = indM(1:num); 
    M(j) = median(ENERGY{j}(lowIND1)); 
             
end
[~,ind] = sort(M); 
% reorder colours based on new plotting order
ModelColors = ModelColors(ind,:); 

I = cell(size(ENERGY,2),1); 

% select lowest values
for j=1:size(ENERGY,2)
    selM = ind(j); 
    switch selWhat
        case 'lowest'
            [v,indM] = sort(ENERGY{selM}); 
        case 'highest'
            [v,indM] = sort(ENERGY{selM}, 'descend'); 
        case 'all'
            % just plot all submitted values - selected before
            indM = 1:length(ENERGY{selM}); 
            % make this variable equal to all existing values
            num = length(ENERGY{selM});
    end

    lowIND2 = indM(1:num); 
    if contains(mtype{selM}, '-')
        mtypeNEW{j} = erase(mtype{selM},'-'); 
    else
        mtypeNEW{j} = mtype{selM}; 
    end
    S.(mtypeNEW{j}) = ENERGY{selM}(lowIND2);
    I{selM} = lowIND2; 
    % the order of I is the same as input - not ordered by energy
end


% put data into structure
if doPlot
fig = figure('color','white');
set(gcf, 'Position', [10 10 1200 500]);

set(gca,'FontSize',18)
subplot(12,4,1:24);
violins = violinplot(S);
ylabel('Model fit (max KS)')

% change colour
violinLabels = cell(size(mtypeNEW,1),1); 

for i=1:size(mtypeNEW,2)
    violins(1,i).ViolinColor = ModelColors(i,:);
    violins(1,i).EdgeColor = ModelColors(i,:);
    violins(1,i).BoxColor = [.25 .25 .25];
   
    violinLabels{i} = mtype{ind(i)}; 
end

xticklabels(violinLabels); 
xtickangle(30)
end

end

