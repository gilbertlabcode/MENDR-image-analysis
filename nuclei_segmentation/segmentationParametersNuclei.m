function p = segmentationParametersNuclei()

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function returns a structure p with parameters for nuclei
% segmentation. High-level parameter tuning should be done with this file
%
% Author: Jose L. Cadavid, University of Toronto, 2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Volume parameters

% All volumes in cubic microns

% Minimum volume allowable for blob
p.vMin = 245;
% Median volume of single nuclei
p.vMedian = 528;
% Rough maximum volume allowable per single nuclei
p.vMax = 1228;
% Number of white neighbours below which a voxel is considered a boundary
% voxel (20-24 is good)
p.nNeigh=20;

%% Penalties for the objective function to determine the optimal number of ellipsoids

% Penalty for image mismatch
p.pM = 0.1;
% Penalty for ellipsoids being outside of size range
p.pV = 0.1;
% Penalty for ellipsoids overlapping
p.pO = 0.05;
% Penalty for aspect ratio (prevent super ellongated ellipsoids)
p.AR = 0.02;

% Parameters calculating penalties based on logistic functions in the form
% [k xm]. xm is the value at which the penalty is 50% and xmin is the value
% at which the penalty is 5%. Based on these two values, k for the logistic
% function is calculated

% For overlap
xm = 20;
xmin = 10;
p.logOver = [-log(1/0.05-1)/(xmin-xm), xm];
% For ellipsoid volume: Lower range - defined in proportion to the minimum
% nuclear size
xm = 0.9*p.vMin;
xmin = 0.7*p.vMin;
p.logVol(1,:) = [-log(1/0.05-1)/(xmin-xm), xm];
% For ellipsoid volume: Upper range - defined in proportion to max nuclear
% size
xm = 1.6*p.vMax;
xmin = 1.0*p.vMax;
p.logVol(2,:) = [-log(1/0.05-1)/(xmin-xm), xm];
% For aspect ratios
xm = 3;
xmin = 2;
p.logAR = [-log(1/0.05-1)/(xmin-xm), xm];
