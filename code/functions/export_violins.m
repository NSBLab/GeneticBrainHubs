function S_exp = export_violins(dataCell)

for rr=1:length(dataCell)
    
    num_points = cellfun(@length,dataCell);
    num_violins = size(dataCell,1);
    S_exp = nan(max(num_points),num_violins);
    
    for pp=1:num_violins
        
        S_exp(1:num_points(pp),pp) = dataCell{pp};
        
    end

end

end