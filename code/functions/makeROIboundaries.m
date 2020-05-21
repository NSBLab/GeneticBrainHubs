
UniqueROIs = unique(cdata_continuous);
UniqueROIs(UniqueROIs==0) = []; 
for i = 1:length(UniqueROIs)
ROIid = UniqueROIs(i);
data = cdata_continuous;
%ROIid = 34;
points=find(data==ROIid);

faces_IND = find(sum(ismember(double(in.faces),points),2) > 0);

ROIfaces = double(in.faces(faces_IND,:));

U = unique(ROIfaces(:));

ROIverts = lh_pialverts(U,:);

ROIfaces = changem(ROIfaces,1:length(U),U);

TR = triangulation(ROIfaces,ROIverts);
F{i} = freeBoundary(TR);

end
