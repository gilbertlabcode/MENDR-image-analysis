function [ImL,idxClust]=labelSplitNuclei(sIm,X,Me,oe)

%Finally cluster the voxels in a blob contained in Im according to the k
%ellipsoids defined by the mappings Me and offsets oe based on minimum distance. Returns a
%labeled image

%Initialize label image
ImL=zeros(sIm);

%Get the distance from every point to the ellipse centers

%Get number of clusters
[k,~]=size(oe);
%Initialize distance matrix
D=zeros(length(X),k);

for i=1:k
    %Ellipsoid radius
    R=(det(Me(:,:,i)))^(1/3);
    %Ellipsoid covariance matrix based on eq (x-xc)'*A(x-xc)
    A=(inv(Me(:,:,i)))^2;
    %Substract center offset
    Xd=X-oe(i,:);
    %Distances (Normalized Mahalanobis) - vectorized to work quickly,
    %should use A' but A is symmetric anyways
    D(:,i)=R*sqrt(sum(Xd.*(Xd*A),2));
%     D(:,i)=(R^2-R^2*(sum(Xd.*(Xd*A),2))).^2;
end

%Find closest cluster based on M distance to center
[~,idxClust]=min(D,[],2);

%Assign labels to voxels in label image
ImL(sub2ind(sIm,X(:,1),X(:,2),X(:,3)))=idxClust;