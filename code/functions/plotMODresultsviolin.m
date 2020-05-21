% make violin plot for modelling results
% make a violin plot
function [S, fig] = plotMODresultsviolin(MEASURE, mtype, measureName)

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
%ModelColors(14,:) = []; 

% order models based on mean fit values
for j=1:size(MEASURE,1)
    M(j) = mean(MEASURE{j}); 
end
[~,ind] = sort(M); 
% reorder colours based on new plotting order
ModelColors = ModelColors(ind,:); 

% select lowest values
for j=1:size(MEASURE,1)
    selM = ind(j); 
    %[v,indM] = sort(MEASURE{selM}); 
    %lowIND2 = indM(1:num); 
    if contains(mtype{selM}, '-')
        mtypeNEW{j} = erase(mtype{selM},'-'); 
    else
        mtypeNEW{j} = mtype{selM}; 
    end

    S.(mtypeNEW{j}) = MEASURE{selM}; 
end


% put data into structure
fig = figure('color','white');
set(gcf, 'Position', [10 10 1200 500]);

set(gca,'FontSize',18)
violins = violinplot(S);
ylabel(sprintf('%s',measureName))

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

