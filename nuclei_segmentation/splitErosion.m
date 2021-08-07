function im2=splitErosion(im,S,Sc)

%%% Erode a blob in im with structuring element S, relabel and dilate EACH
%%% blob (that would be the difference between this and opening)

[M,N,P]=size(im);
%Erode image to see if it is split
im=imerode(im,S);
%Relabel new image
L=bwlabeln(im);
%Get voxels in each region
stats=regionprops3(L,'VoxelIdxList');

%Dilate each new blob if there is more than one. Some holes might have been
%opened so we perform closing in each blob (with a smaller structuring
%element, if possible, to preserve boundary infoin blobs)

%Output image with new labels
im2=zeros(size(im));

%Dilate each image (to process faster take bounding box)
for i=1:height(stats)
    
    %Get bounding box of blob
    [x,y,z]=ind2sub(size(im),stats.VoxelIdxList{i});
    %Create extra padding of 4 (should be a function of the structuring
    %element, but I am using radius 3) because blob will be dilated
    xmin=max(1,min(x)-4);
    xmax=min(M,max(x)+4);
    ymin=max(1,min(y)-4);
    ymax=min(N,max(y)+4);
    zmin=max(1,min(z)-4);
    zmax=min(P,max(z)+4);
    
    %Create dummy image that is bounding box of blob
    imDummy=zeros(size(im));
    imDummy(stats.VoxelIdxList{i})=1;
    %Dilate blob
    imDummy=imdilate(imDummy(xmin:xmax,ymin:ymax,zmin:zmax),S);
    
    %Close blob
%     if imEuler3d(imDummy,6)<1
%         imDummy=imclose(imDummy,Sc);
%     end
    %Put this info in the big image
    %Find coordinates where voxels are non-zero with respect to to bounding box
    [xp, yp, zp]=ind2sub([xmax-xmin+1,ymax-ymin+1,zmax-zmin+1],find(imDummy(:)==1));
    %Translate these coordinates to the coordinate system of the box. Add
    %the new label. Add 1 to these voxels. That way, things that were
    %already labelled (overlaps) will have a value > 1
    im2(sub2ind(size(im),xp+xmin-1,yp+ymin-1,zp+zmin-1))=i; %im2(sub2ind(size(im),xp+xmin-1,yp+ymin-1,zp+zmin-1))+
    
end