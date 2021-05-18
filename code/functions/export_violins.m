function S_exp = export_violins(dataCell)

for rr=1:length(dataCell)
    
    num_points = cellfun(@length,dataCell(rr));
    num_violins = size(dataCell{rr},1);
    d_exp = dataCell{rr};
    S_exp = nan(max(num_points),num_violins);
    
    for pp=1:num_violins
        
        S_exp(1:num_points(pp),pp) = d_exp(pp);
        
    end

end

end