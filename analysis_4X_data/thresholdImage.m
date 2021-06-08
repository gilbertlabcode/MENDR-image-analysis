function [pArea, tOpt, flag, param] = thresholdImage(img, plt)

% thresholdImage takes in a grayscale image (8-bit) and attempts to find an optimal
% threshold to binarize it by assuming that its intensity histogram is a combination of
% two log-normal distributions (foreground and background). The histogram is decomposed
% with a least-squares method, and the threshold is defined as the intensity at which
% class uncertainty is maximum.
%
% USAGE:
%
%   [pArea, tOpt, flag, param] = thresholdImage(img, plt)
%
% INPUTS:
%
%   img:    8-bit grayscale image to be thresholded
%   plt:    Optional boolean that indicates whether intermediate plots should be
%           made (1) or not (0). Default is 0
%
% OUTPUTS:
%
%   pArea:  Percentage coverage of foreground pixels
%   tOpt:   Optimal intensity threshold
%   flag:   Boolean that indicates whether the two classes of pixels are
%           not properly separated (flag = 1; threshold not reliable)
%   param:  Parameters of the fitted log-normal distributions
%
% AUTHOR:
%   Jose L. Cadavid, University of Toronto, 2021

%% Get histogram (CDF) of image and define bin centers for passing to the parameter fit routine

% Set plot to 0 by default if argument is not included
if nargin == 1
    plt = 0;
end

%Get counts and edges of bins
[counts, edges] = histcounts(img,'normalization','cdf');
%Get centers of bins
centers = 0.5*(edges(2:end)+edges(1:end-1));

%% Fit CDF to two log normal distirbutions and get optimal threshold tOpt
[param, tOpt, flag] = fitCDF(centers, counts, plt);

%% Get percentage coverage based on optimal threshold
pArea = 100*sum(img(:)>=tOpt)/numel(img);
