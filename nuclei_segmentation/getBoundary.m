function imgBound = getBoundary(img, nNeigh)

% getBoundary gets the boundary voxels of blobs in a binary image with many blobs.
% a boundary voxel is defines as a white voxel with fewer than 24
% neighbours (the number can be tuned, maximum 26). Number of neghbours
% counted quickly with convolution. Return an image of those boundaries.
%
% USAGE:
%
%   imgBound = getBoundary(img, nNeigh)
%
% INPUTS:
%
%   img:       Binary 3D image
%   nNeigh:    Number of voxels required in the neighborhood
%
% OUTPUTS:
%
%   imgBound:    Image of boundaries
%
% AUTHOR:
%   Jose L. Cadavid, University of Toronto, 2021

% Convolution kernel for counting neighbours
F = ones(3,3,3);
F(2,2,2) = 0;
% Convolve: Voxels in new image have a value equal to the number of white
% neighbours
imgNeigh = convn(double(img),F,'same');
% Eliminate voxels from edges of the images for calculating the boundary
% voxels of the blob
img(1,:,:) = 0;
img(end,:,:) = 0;
img(:,:,end) = 0;
img(:,:,1) = 0;
img(:,1,:) = 0;
img(:,end,:) = 0;

% Classify a voxel as a boundary voxel depending on it being white and
% having less than nNeigh white neighbours (can be changed, but 20-24 is a good
% number)
imgBound = (img > 0) & (imgNeigh < nNeigh);