function makeROIcont(in)

NonContinuousVertices = [];
cdata_continuous = in.cdata;
parcdata = cdata_continuous;

UniqueROIs = unique(in.cdata);

for i = UniqueROIs

ROIid = i;
points=find(parcdata==ROIid);

faces_IND = find(sum(ismember(in.faces,points),2) > 0);

ROIfaces = lh_faces(faces_IND,:);

ROIfaces = ROIfaces.*ismember(ROIfaces,points);

U = unique(ROIfaces(:));

ROIfaces = changem(ROIfaces,1:length(U),U);

        A = zeros(length(U));
   for j = 1:length(faces_IND)

       X = ROIfaces(j,1);
       Y = ROIfaces(j,2);
       Z = ROIfaces(j,3);

       A(X,Y) = 1;
       A(Y,X) = 1;
       A(Z,Y) = 1;
       A(Y,Z) = 1;
       A(X,Z) = 1;
       A(Z,X) = 1;

   end
      S = sum(A);
      A(S<2,:) = 0;
      A(:,S<2) = 0;

      A(1,:) = [];
      A(:,1) = [];
      U(1) = [];
   comps = graphComponents(A);

   NonContinuous = find(comps ~= mode(comps));

   NonContinuousVertices = [NonContinuousVertices U(NonContinuous)'];

end

while ~isempty(NonContinuousVertices)

    for j = 1:length(NonContinuousVertices)
        VertID = NonContinuousVertices(j);
        currentROIid = lh_aparc(VertID);
        faces_IND = find(sum(ismember(lh_faces,VertID),2) > 0);
        faces2check = lh_faces(faces_IND,:);
        neiborVertxROIids = lh_aparc(unique(faces2check(:)));
        neiborVertxROIids(neiborVertxROIids==currentROIid) = NaN;
        LargestNeighbor = mode(neiborVertxROIids);
        %if ~isnan(LargestNeighbor)
            cdata_continuous(VertID) = LargestNeighbor;
            NonContinuousVertices(j) = NaN;
        %end
    end

    NonContinuousVertices(isnan(NonContinuousVertices)) = [];

end
end
