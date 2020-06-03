% make violin plot for modelling results
% make a violin plot
function [S,INDselected] = plotMODviolin(ENERGY, num, mtype, selWhat, doPlot)
if nargin<5
    doPlot = true; 
end

S = struct; 
% assign colours based on original order

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
    INDselected{selM} = indM(1:num); 
    if contains(mtype{selM}, '-')
        mtypeNEW{j} = erase(mtype{selM},'-'); 
    else
        mtypeNEW{j} = mtype{selM}; 
    end
    S.(mtypeNEW{j}) = ENERGY{selM}(lowIND2); 
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

