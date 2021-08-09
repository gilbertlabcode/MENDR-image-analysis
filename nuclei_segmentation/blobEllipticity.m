function E = blobEllipticity(X, Vol, Area)

% blobEllipticity calculates a shape metric similar to sphericity that
% we call "ellipticity". This metric is a generalization of sphericity
% that uses the moments of inertia of a blob to approximate it as an
% ellipsoid.
%
% USAGE:
%
%   E = blobEllipticity(X, Vol, Area)
%
% INPUTS:
%
%   X:      Matrix (nx3) of coordinates of voxels that makeup blob
%   Vol:    Volume of blob (in voxels)
%   Area:   Surface area of blob
%
% OUTPUTS:
%
%   E:      Ellipticity of blob (1 is perfect ellipsoid)
%
% AUTHOR:
%   Jose L. Cadavid, University of Toronto, 2020

% Get covariance matrix from which the ellipsoid that describes the data is
% obtained 
covMat=cov(X);
% Get eigenvalues of covariance matrix, which are related to the semiaxis of
% the ellipsoid
[~,d] = eig(covMat);
% Get the radii of the ellipse sorted in descending order. Formula from
% https://imagejdocu.tudor.lu/tutorial/plugins/3d_ellipsoid.
R = sort(sqrt(5*diag(d)),'descend');
% Normalize the radii so that the volume of the ellipse matches the volume
% of the blob (scaling of all radii)
R = R*(Vol/(4*pi*prod(R)/3))^(1/3);
% Get the surface area of an ellipsoid with the computed radii (the fitted
% ellipsoid has the same volume as the 3D blob by definition). Using the
% elliptical integrals
theta = acos(R(3)/R(1));
k = sqrt((R(1)^2*(R(2)^2-R(3)^2)/(R(2)^2*(R(1)^2-R(3)^2))));
Ae = 2*pi*R(3)^2+2*pi*R(1)*R(2)/sin(theta)*(ellipticE(theta,k)*(sin(theta))^2+ellipticF(theta,k)*(cos(theta))^2);
% Calculate ellipticity E as the ratio of the surface area of the fitted
% ellipsoid vs the surface are of the blob. Sphericity drops below 1 for
% flattened ellipsoids, but this metric does not
E = Ae/Area;


