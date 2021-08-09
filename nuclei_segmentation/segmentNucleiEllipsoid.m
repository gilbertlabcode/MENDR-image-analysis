function [imgLabelNew,dNuc,oeNuc] = segmentNucleiEllipsoid(imgNuc, scale)

% segmentNucleiEllipsoids takes a 3D image (isotropic) of nuclei and
% attempts to segment it. Segmentation on fitting ellipsoids to
% binary blobs and k-means with Mahalanobis distance. Initial binarization 
% done by smoothing the image with a non-linear filter to preserve edges as 
% much as possible, and then aplying Otsu's method
%
% USAGE:
%
%   [imgLabelNew,dNuc,oeNuc] = segmentNucleiEllipsoid(imgNuc, scale)
%
% INPUTS:
%
%   imgNuc:    3D double grayscale image of the nuclei
%   scale:     Size of a voxel in microns
%
% OUTPUTS:
%
%   nucStats:    Table with volume of segmented nuclei. It also indicates
%                if nuclei are GFP, SAA or double-positive
%
% AUTHOR:
%   Jose L. Cadavid, University of Toronto, 2021

%% Get segmentation parameters
PARAM = segmentationParametersNuclei();

%% Step 1: Binarize image by using global threshold and clean up for segmentation

% Blur
imgFilt = imdiffusefilt(imgNuc, 'NumberOfIterations', 3, 'GradientThreshold', 0.2*max(imgNuc(:)));
% Binarize
imgBin = imbinarize(imgFilt/255);

% Label connected components and remove small blobs (of less than Vmin)
% Note: Tuning point, volume of blobs to be removed can be changed

% Get volume of blobs
stats = regionprops3(bwlabeln(imgBin), 'Volume', 'VoxelIdxList');
% Transform that volume from voxels to cubic microns with the scale provided
stats.Volume = stats.Volume*scale^3;
for i = 1:height(stats)
    if stats.Volume(i) < PARAM.vMin
        % Small component, delete voxels
        imgBin(stats.VoxelIdxList{i,1}) = 0;
    end
end

% Apply erosion and dilation to split blobs where there are small
% connecting regions. This also partially gets rid of tiny objects. By
% splitting large blobs with spherical structuring elements we can reduce
% computation time in clustering; especially if we split huge blobs.
% Small holes are also filed with morphological closing if necessary

% Note: The size of the structuring elements can be tweaked or made to be
% adaptive
imgBin = splitErosion(imgBin, strel('sphere',3));
% Optional: Remove objects that touch border of image. For analysing the
% morphology of nuclei, these objects should be ignored
imgBin = imclearborder(imgBin>0).*imgBin;

% Note: At this point, the binary image is not longer binary, but the different
% split blobs are labelled. We go through a second round of removing tiny
% objects that might have appeared. We also keep a list of the blobs that
% are removed for relabelling

% Get volume of blobs
stats = regionprops3(imgBin,'Volume','VoxelIdxList');
% Transform that volume from voxels to cubic microns with the scale provided
stats.Volume = stats.Volume*scale^3;
nBlobs = height(stats);

% Tag blobs for eliminating
list = stats.Volume < PARAM.vMin;
imgLabel = imgBin; 
% Relabel without the indices of the removed blobs
ii=1;
for i = 1:nBlobs
    if list(i) % eliminate
        imgLabel(stats.VoxelIdxList{i,1}) = 0;
    else % relabel
        imgLabel(stats.VoxelIdxList{i,1}) = ii;
        %Next label
        ii = ii+1;
    end
end

%% Step 2: Split touching nuclei in the blobs with a modified k-means algorithm with Mahalanobis distance

% Relabel image: At this point some adjacent nuclei have been labelled with
% different numbers in the image erosion algorithm, but if we relabel them
% they will be fused because they are connected (the erosion algorithm
% doesn't always create a black space between them because the dilation can
% make them touch again, which is why opening doesn't work the same way).

% Get relevant properties of blobs 
statsL = regionprops3(imgLabel,'Volume','VoxelIdxList','SurfaceArea');
% Scale volume from voxels to microns
statsL.Volume = statsL.Volume*scale^3;
%Scale area from pix square to micron 2
statsL.SurfaceArea = statsL.SurfaceArea*scale^2;
nBlobs = height(statsL);

% Get the boundary voxels of the blobs by convolution and counting
% neighbours. This function returns it with no labels (only 0 and 1), so we
% multiply it with the fully labelled image to label those boundary voxels
imgBound = getBoundary(imgLabel>0,PARAM.nNeigh).*imgLabel;

%Initialize a new label image
imgLabelNew = zeros(size(imgLabel));

% Matrices for storing ellipsoid centers and semiaxes 
dNuc = [];
oeNuc = [];
%Counter for labelled nuclei
nucCount = 0;

% Split blobs, one by one, based on ellipsoid clustering
for i = 1:nBlobs
    
    % Get optimal split of blob in terms of ellipsoids (optimal number
    % calculated internally). Also return the center (oeOpt) and semiaxes
    % length (dOpt) of those ellipsoids. Ignoring the rotation angles but
    % could be obtained as well
    [dDummy,oeDummy,idxVox,idxClust] = splitNucleiClusteringFaster(imgLabel==i,...
                                        imgBound==i, statsL(i,:), PARAM, scale);
    % Get the number of new ellipsoids
    [kOpt,~] = size(oeDummy);
    % Label the voxels of the blob indicated by idxVox with the labels
    % idxClust (adding the current nuclei counter)
    imgLabelNew(idxVox) = idxClust + nucCount;
    % Store semiaxes length and centers
    dNuc = [dNuc;dDummy];
    oeNuc = [oeNuc;oeDummy];
    % Update nuclei counter
    nucCount = nucCount + kOpt;
end

