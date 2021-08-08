function imBound=getBoundary(Im,nNeigh)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Get the boundary voxels of blobs in a binary image with many blobs.
%%% a boundary voxel is defines as a white voxel with fewer than 24
%%% neighbours (the number can be tuned, maximum 26). Number of neghbours
%%% counted quickly with convolution. Return an image of those boundaries
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Convolution kernel for counting neighbours
F=ones(3,3,3);
F(2,2,2)=0;
%Convolve: Voxels in new image have a value equal to the number of white
%neighbours
imNeigh=convn(double(Im),F,'same');
%Eliminate voxels from edges of the images for calculating the boundary
%voxels of the blob
Im(1,:,:)=0;
Im(end,:,:)=0;
Im(:,:,end)=0;
Im(:,:,1)=0;
Im(:,1,:)=0;
Im(:,end,:)=0;

%Classify a voxel as a boundary voxel depending on it being white and
%having less than nNeigh white neighbours (can be changed, but 20-24 is a good
%number)
imBound=(Im>0) & (imNeigh<nNeigh);