function imgSplit = splitErosion(img,S)

% splitErosion takes a binary image and attempts to break blobs by
% performing erosion, relabelling, and dilation of new labels (the
% relabelling makes this different from morphological opening).
%
% USAGE:
%
%   imgSplit = splitErosion(img,S)
%
% INPUTS:
%
%   img:    Binary image
%   S:      Spherical structuring element created with strel()
%
% OUTPUTS:
%
%   imgSplit:    Label image of split blobs
%
% AUTHOR:
%   Jose L. Cadavid, University of Toronto, 2021

% Get image size
[M,N,P]=size(img);
%% Erode image 
img = imerode(img,S);
% Relabel new image
imgLab = bwlabeln(img);
% Get voxels in each labelled blob
stats = regionprops3(imgLab,'VoxelIdxList');
nBlobs = height(stats);

%% Dilate each new blob if there is more than one

% Output image with new labels
imgSplit=zeros(size(img));

% Dilate each each blob. Make faster by only taking its bounding box
for i = 1:nBlobs
    %Get bounding box of blob
    % Get x y z coordinates of voxels in blob
    [x,y,z] = ind2sub(size(img), stats.VoxelIdxList{i});
    % Create extra padding of 4 (should be a function of the structuring
    % element, but I am using radius 3) because each blob will be dilated
    xmin = max(1,min(x)-4);
    xmax = min(M,max(x)+4);
    ymin = max(1,min(y)-4);
    ymax = min(N,max(y)+4);
    zmin = max(1,min(z)-4);
    zmax = min(P,max(z)+4);
    
    % Create dummy image that is bounding box of blob
    imgDummy = zeros(size(img));
    imgDummy(stats.VoxelIdxList{i}) = 1;
    % Dilate blob
    imgDummy = imdilate(imgDummy(xmin:xmax,ymin:ymax,zmin:zmax),S);
    
    % Put this info in the big image
    % Find coordinates where voxels are non-zero with respect to to bounding box
    [xp, yp, zp] = ind2sub([xmax-xmin+1,ymax-ymin+1,zmax-zmin+1],find(imgDummy(:)==1));
    % Translate these coordinates to the coordinate system of the box. Add
    % the new label. Add 1 to these voxels. That way, things that were
    % already labelled (overlaps) will have a value > 1
    imgSplit(sub2ind(size(img),xp+xmin-1,yp+ymin-1,zp+zmin-1)) = i;     
end