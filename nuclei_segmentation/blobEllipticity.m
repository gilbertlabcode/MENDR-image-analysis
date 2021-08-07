function E=blobEllipticity(X,Vol,Area)

%%% Written by Jose Cadavid, University of Toronto, 2019

%%% Get the "ellipticity" of the 3D object in the binary image Mask. Using
%%% voxel locations given by matrix X; the blob has volume Vol and surface
%%% area Area.

%%% This metric is a generalization of sphericity for objects that resemble
%%% ellipsoids, and I don't know if it exists already and I don't know it's
%%% name


%Get covariance matrix from which the ellipsoid that describes the data is
%obtained (different from algebraic fitting an ellipsoid, this is faster
%for a single ellipsoid)
Mat=cov(X);
%Get eigenvalues of covariance matrix, which are related to the semiaxis of
%the ellipsoid
[~,d]=eig(Mat);
%Get the radii of the ellipse sorted in descending order. I don't know why there is a 5 but this
%equation works https://imagejdocu.tudor.lu/tutorial/plugins/3d_ellipsoid.
%I suspect it comes from the moment of inertia tensor of the ellipsoid
%(second moment), which has a 1/5 in there https://www.wolframalpha.com/input/?i=ellipsoid
R=sort(sqrt(5*diag(d)),'descend');
%Normalize the radii so that the volume of the ellipse matches the volume
%of the blob (scaling of all radii)
R=R*(Vol/(4*pi*prod(R)/3))^(1/3);
%Get the surface area of an ellipsoid with the computed radii (the fitted
%ellipsoid has the same volume as the 3D blob by definition). Using the
%elliptical integrals
theta=acos(R(3)/R(1));
k=sqrt((R(1)^2*(R(2)^2-R(3)^2)/(R(2)^2*(R(1)^2-R(3)^2))));
Ae=2*pi*R(3)^2+2*pi*R(1)*R(2)/sin(theta)*(ellipticE(theta,k)*(sin(theta))^2+ellipticF(theta,k)*(cos(theta))^2);
%Calculate ellipticity E as the ratio of the surface area of the fitted
%ellipsoid vs the surface are of the blob. Sphericity drops below 1 for
%flattened ellipsoids, but this metric does not
E=Ae/Area;


