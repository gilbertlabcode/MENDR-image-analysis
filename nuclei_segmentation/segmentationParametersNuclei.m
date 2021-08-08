function PARAM = segmentationParametersNuclei()

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function returns a structure PARAM with parameters for nuclei
% segmentation. High-level parameter tuning should be done with this file
%
% Author: Jose L. Cadavid, University of Toronto, 2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Volume parameters

% All volumes in cubic microns

% Minimum volume allowable for blob
PARAM.vMin = 245;
% Median volume of single nuclei
PARAM.vMedian = 528;
% Rough maximum volume allowable per single nuclei
PARAM.vMax = 1228;
% Number of white neighbours below which a voxel is considered a boundary
% voxel (20-24 is good)
PARAM.nNeigh=20;

%% Penalties for the objective function to determine the optimal number of ellipsoids

% Penalty for image mismatch
PARAM.pM = 0.1;
% Penalty for ellipsoids being outside of size range
PARAM.pV = 0.1;
% Penalty for ellipsoids overlapping
PARAM.pO = 0.05;
% Penalty for aspect ratio (prevent super ellongated ellipsoids)
PARAM.AR = 0.02;

% Parameters calculating penalties based on logistic functions in the form
% [k xm]. xm is the value at which the penalty is 50% and xmin is the value
% at which the penalty is 5%. Based on these two values, k for the logistic
% function is calculated

% For overlap
xm = 20;
xmin = 10;
PARAM.logOver = [-log(1/0.05-1)/(xmin-xm), xm];
% For ellipsoid volume: Lower range - defined in proportion to the minimum
% nuclear size
xm = 0.9*PARAM.vMin;
xmin = 0.7*PARAM.vMin;
PARAM.logVol(1,:) = [-log(1/0.05-1)/(xmin-xm), xm];
% For ellipsoid volume: Upper range - defined in proportion to max nuclear
% size
xm = 1.6*PARAM.vMax;
xmin = 1.0*PARAM.vMax;
PARAM.logVol(2,:) = [-log(1/0.05-1)/(xmin-xm), xm];
% For aspect ratios
xm = 3;
xmin = 2;
PARAM.logAR = [-log(1/0.05-1)/(xmin-xm), xm];
