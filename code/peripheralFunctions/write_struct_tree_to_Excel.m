function write_struct_tree_to_Excel(xlsfile,strct)
% function write_struct_tree_to_Excel(xlsfile)
%   inputs: name of Excel file
%           struct (which might have other structs)
%   output: none

%% Write Struct Data
txt = fieldnames(strct) ;
sel = ones(size(txt)) ;
for i = 1:length(txt)
    sel(i) = isstruct(strct.(txt{i})) ;
end

i_not_struct = find(~sel) ;
i_struct = find(sel) ;

x = [fieldnames(strct) struct2cell(strct)] ;
xlswrite(xlsfile ,x(i_not_struct,:),1,'a1') ; % winopen(xlsfile)

loc1 = size(x(i_not_struct,:),1) + 1 ;
loc_txt = {sprintf('a%d',loc1),sprintf('b%d',loc1)} ;
for i = 1:length(i_struct)
    xlswrite(xlsfile ,x(i_struct(i),1),1,loc_txt{1}) ; % output name of the struct
    sub_strct = x{i_struct(i),2} ; % get actual struct
    y = [fieldnames(sub_strct) struct2cell(sub_strct)] ;
    for j = 1:size(y,1)
        is_array = numel(y{j,2}) > 1 ;
        if is_array
            txt = sprintf('%g,',y{j,2}) ; txt(end) = [] ;
            y{j,2} = txt ;
        end
    end
    xlswrite(xlsfile ,y,1,loc_txt{2}) ;
    loc1 = loc1 + size(y,1) ;
    loc_txt = {sprintf('a%d',loc1),sprintf('b%d',loc1)} ;
end
