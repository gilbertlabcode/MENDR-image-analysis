function [nucStats, nNuclei, nGFP, nSAA, nGFPSAA] = getNuclearFractions(imgNuc, imgGFP, imgSAA, scale)

% getNuclearFractions takes three isotropic 3D grayscale images (nuclei, GFP, SAA) and
% calculates whether nuclei (after segmentation) are within a GFP or SAA
% positive structure
%
% USAGE:
%
%   [nucStats, nNuclei, nGFP, nSAA, nGFPSAA] = getNuclearFractions(imgNuc, imgGFP, imgSAA, scale)
%
% INPUTS:
%
%   imgNuc:    3D 8-bit grayscale image of the nuclei
%   imgGFP:    3D 8-bit grayscale image of the GFP channel
%   imgEdU:    3D 8-bit grayscale image of the EdU channel
%   scale:     Size of a voxel in microns
%
% OUTPUTS:
%
%   nucStats:    Table with volume of automatically segmented nuclei. It also indicates
%                if nuclei are GFP, SAA or double-positive
%   nNuclei:     Number of detected nuclei
%   nGFP:        Number of detected nuclei in a GFP structure
%   nSAA:        Number of detected nuclei in a SAA structure
%   nGFPSAA:     Number of detected nuclei in a SAA+GFP+ structure
%
% AUTHOR:
%   Jose L. Cadavid, University of Toronto, 2021

% Convert nuclei image to double, and GFP and SAA to unit8
imgNuc = double(imgNuc);
imgGFP = uint8(imgGFP);
imgSAA = uint8(imgSAA);

%% Binarize and segment nuclei with the iterative k-cluster
[imgLabel, ~, ~, imgBig] = segmentNucleiEllipsoid(imgNuc, scale);
% Get properties of nuclei (only volume)
nucStats = regionprops3(imgLabel, 'Volume', 'VoxelIdxList');
% Remove any nuclei with zero volume that might have been mislabeled
nucStats(nucStats.Volume(:)==0,:) = [];
nNuclei = height(nucStats);
meanVol = mean(nucStats.Volume(:));

%% Binarize GFP and SAA channels

% Heavily filter GFP and SAA with a median filter
imgGFP = medfilt3(imgGFP,[5 5 5]);
imgSAA = medfilt3(imgSAA,[5 5 5]);

% Threshold these images using a MATLAB implementation of the triangle
% method
% (https://www.mathworks.com/matlabcentral/fileexchange/28047-gray-image-thresholding-using-the-triangle-method)

%Get histograms to feed the triangle method
[GFPhist, ~] = imhist(imgGFP);
[SAAhist, ~] = imhist(imgSAA);

%Get threshold
tGFP = triangle_th(GFPhist,256)*max(imgGFP(:));
tSAA = triangle_th(SAAhist,256)*max(imgSAA(:));

% To define whether a nucleus is inside the muscle fiber we say that at least
% 75% of its volume must be masked in the GFP or SAA channels (i.e. 75%
% of the pixels in the nucleus must be brighter than the thresholds).
% This is to allow for small holes that appear during threshold or
% imperfect color alignment from the confocal and can be tuned to be more
% or less stringent

for i = 1:nNuclei
    nucRegion = nucStats.VoxelIdxList{i,1};
    nucVol = nucStats.Volume(i);
    % Number of pixels in the region brighter than the threshold
    nGFP = sum(imgGFP(nucRegion) > tGFP);
    % Number of pixels in the region brighter than the threshold
    nSAA = sum(imgSAA(nucRegion) > tSAA);
    % Is the nucleus considered GFP or SAA?
    % Add a column to the NucStats table with 1 or 0
    nucStats.GFP(i) = (nGFP/nucVol >= 0.75);
    nucStats.SAA(i) = (nSAA/nucVol >= 0.75);
    nucStats.SAAGFP(i) = (nGFP/nucVol >= 0.75) && (nSAA/nucVol >= 0.75);
end

% Get number of nuclei, number of GFP, SAA and doublepos nuclei
nGFP = sum(nucStats.GFP(:));
nSAA = sum(nucStats.SAA(:));
nGFPSAA = sum(nucStats.SAAGFP(:));

%% Process the big blobs

% These blobs are not segmented properly so we estimate the number of
% nuclei the have by dividing their volume by the mean volume of segmented
% nuclei. We then estimate the number of GFP/SAA based on the fraction of
% SAA/GFP voxels times the number of estimated nuclei. In practice, few
% blobs will go through this process so the impact of this estimation step
% is not huge

bigStats = regionprops3(imgBig, 'Volume', 'VoxelIdxList');
for i = 1:height(bigStats)
    nucRegion = bigStats.VoxelIdxList{i,1};
    nucVol = bigStats.Volume(i);
    % Estimate number of nuclei based on volume
    nNew = round(nucVol/meanVol);
    % Estimate the fraction of the blob that overlaps with GFP/SAA
    fracGFP = sum(imgGFP(nucRegion) > tGFP)/nucVol;
    fracSAA = sum(imgSAA(nucRegion) > tSAA)/nucVol;
    fracSAAGFP = sum(imgGFP(nucRegion) > tGFP & imgSAA(nucRegion) > tSAA)/nucVol;
    % Add counts to total counts 
    nNuclei = nNuclei + nNew;
    nGFP = nGFP + round(fracGFP*nNew);
    nSAA = nSAA + round(fracSAA*nNew);
    nGFPSAA = nGFPSAA + round(fracSAAGFP*nNew);
end
    



