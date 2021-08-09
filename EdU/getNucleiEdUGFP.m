function [nucGFP, nGFP, nEdUGFP] = getNucleiEdUGFP(imgNuc, imgGFP, imgEdU, scale, vMin)

% getNucleiEdUGFP takes three isotropic 3D grayscale images (nuclei, gfp, EdU) and
% calculates the number of GFP objects that have EdU+ nuclei. It does so by
% binarizing the images with automatic thresholds.
%
% USAGE:
%
%   [threshGFP, nucGFP, nGFP, nEdUGFP] = getNucleiEdUGFP(imgNuc, imgGFP, imgEdU, scale, vMin)
%
% INPUTS:
%
%   imgNuc:    3D 8-bit grayscale image of the nuclei
%   imgGFP:    3D 8-bit grayscale image of the GFP channel
%   imgEdU:    3D 8-bit grayscale image of the EdU channel
%   scale:     Size of a voxel in microns
%   vMin:      Smallest volume of nuclei in cubic microns
%
% OUTPUTS:
%
%   nucGFP:   Number of nuclei detected in GFP objects
%   nGFP:     Number of GFP objects detected
%   nEdUGFP:  Number of GFP objects with positive EdU signal in them
%
% AUTHOR:
%   Jose L. Cadavid, University of Toronto, 2021

%% Filter and binarize the nuclei channel

% Blur: Median filter
imgNuc = double(medfilt3(imgNuc,[3 3 3]));
% Binarize
imgNucBin = imbinarize(imgNuc/255);

% Label connected components and remove small blobs (of less than Vmin)
% Note: Tuning point, volume of blobs to be removed can be changed

%Get volume of blobs
stats = regionprops3(bwlabeln(imgNucBin),'Volume','VoxelIdxList');
%Transform that volume from voxels to cubic microns with the scale provided
stats.Volume = stats.Volume*scale^3;
for i=1:height(stats)
    if stats.Volume(i) < vMin
        %Small component, delete voxels
        imgNucBin(stats.VoxelIdxList{i,1}) = 0;
    end
end

% At this point we have a clean nuclear mask. We do not attempt to identify
% single nuclei yet.

%% Filter and binarize the GFP channel
imgGFP = double(medfilt3(imgGFP,[3,3,3]));

% Note: Since the GFP fibers are sparse, thresholding methods might identify a lot
% of background. To get around this, I propose using Otsu's method to get an
% initial idea where the fibres are (or any fixed threshold, 20 seems to be
% a good value).
% Otsu's - some images are almost empty so we make sure that the threshold
% is at least a minimum value (determined based on images as 20)
threshMin = 20;
threshGFP = max(threshMin,255*graythresh(imgGFP/255));
imgGFPBin = imgGFP > threshGFP;

% If the threshold was set to less than the minimum it is likely that there are
% no fibers here, so do not rebinarize to avoid error
if threshGFP > threshMin
  % Get distance map: Distances to all black voxels to ANY surface of the GFP
  % blobs
  imgDist = bwdist(imgGFPBin);

  % Get penalties: Based on normal distribution, square difference of GFP
  % intensity vs mean GFP intensity of initial mask, normalized by variance of GFP
  % in mask. This might cut some of the brightest pixels but we will include
  % them again. There is another penalty based on the distance map and an
  % normal distribution normalized by 5 pixel distance (tunable)
  meanGFP = mean(imgGFP(imgGFP(:) > threshGFP));
  varGFP = var(imgGFP(imgGFP(:) > threshGFP));
  score = exp(-(imgGFP - meanGFP).^2/(2*varGFP)).*exp(-imgDist.^2/5^2);

  % Apply threshold based on penalty (the larger the more stringent), 0.2
  % seems fine. Make sure to include again the bright pixels that had passed the
  % initial threshold
  imgGFPBin = (score > 0.2) | imgGFPBin;
end

% Remove GFP objects smaller than about a nucleus
stats = regionprops3(bwlabeln(imgGFPBin),'Volume','VoxelIdxList');
stats.Volume = stats.Volume*scale^3;

% Sanity check, sometimes there is only ONE voxel that passes and then it messes up
% the bracket indexing for the voxel list
if all(stats.Volume <= vMin)
    imgGFPBin = zeros(size(imgGFPBin));
else
    for i = 1:height(stats)
        imgGFPBin(stats.VoxelIdxList{i}) = stats.Volume(i) > vMin;
    end
end

%% Mask the binarized nuclei with the GFP
imgNucGFP = imgGFPBin & imgNucBin;

% Remove any spurious chunk of nuclei that might have been included in the GFP mask
stats = regionprops3(bwlabeln(imgNucGFP),'Volume','VoxelIdxList');
stats.Volume = stats.Volume*scale^3;

%Sanity check, sometimes there is ONE voxel that passes and then it messes up the bracket indexing for the voxel list
if all(stats.Volume <= vMin)
    imgNucGFP = zeros(size(imgNucGFP));
else
    for i = 1:height(stats)
        imgNucGFP(stats.VoxelIdxList{i}) = stats.Volume(i) > vMin;
    end
end

% At this stage overlaps in nuclei are unlikely since they are sparse, so we
% don't split blobs or anything

%% Determine EdU+ nuclei within the GFP objects

%Filter EdU
imgEdU = double(medfilt3(imgEdU,[3 3 3]));

% To determine the nuclei that are EdU+, we get the mean EdU value in the
% region defined by each GFP+ nucleus. Paper fibers are  sometimes visible
% in the EdU channel at a level between background and actual EdU.

% To determine a value to define a nucleus as EdU+ we could use a fixed value
% (say 10) but we can define the background level of EdU by analysing its
% intensity in the region where no nuclei are detected! We can build a
% histogram of EdU pixels in nuclei- regions and define a threshold as the
% value where the cumulative distribution of those pixels is 97.5%

h = imgEdU(~imgNucBin(:));
N = sum(~imgNucBin,'all');
% Get that 97.5% threshold value. This percentile worked well on our images
threshEdU = 0.5;
stop = 0;
while ~stop
    stop = sum(h<threshEdU)/N > 0.975;
    threshEdU = threshEdU + 0.5;
end

%Get mean EdU value in nuclei in the GFP objects
stats = regionprops3(bwlabeln(imgNucGFP),'Volume','VoxelIdxList');
EdUNuc = zeros(height(stats),1);

for i = 1:height(stats)
    EdUNuc(i) = mean(imgEdU(stats.VoxelIdxList{i}));
end

% EdU+ nuclei in the GFP+ objects
nEdUGFP = sum(EdUNuc > threshEdU);
% Total nuclei in the GFP+ objects
nucGFP = height(stats);
% Total number of GFP+ objects
labelGFP = bwlabeln(imgGFPBin);
nGFP = max(labelGFP(:));
