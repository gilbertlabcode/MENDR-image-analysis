function [imgLab, idxClust] = labelSplitNuclei(sizeIm, X, Me, oe)

% labelSplitNuclei clusters the voxels in a blob contained in an image according to the k
% ellipsoids defined by the mappings Me and offsets oe based on minimum distance. Returns a
% labeled image
%
% USAGE:
%
%   [ImL, idxClust] = labelSplitNuclei(sizeIm, X, Me, oe)
%
% INPUTS:
%
%   sizeIm:    Size of the original binary image where the blob is found
%   X:         Matrix with coordinates of the voxels in the blob to be
%              segmented
%   Me:        3x3xk matrix containing shape parameters of each k ellipsoid
%   oe:        3xk matrix containing centroids of eack of k ellipsoids
%
% OUTPUTS:
%
%   imgLab:    Label image where each voxel is labelled based on their
%              proximity to the fitted ellipsoids
%   idxClust:  Vector with indices indicating which cluster each voxel
%              in X belongs to.
%
% AUTHOR:
%   Jose L. Cadavid, University of Toronto, 2021

% Initialize label image
imgLab = zeros(sizeIm);

% Get the distance from every point to the ellipse centers
% Get number of clusters
[k,~] = size(oe);
% Initialize distance matrix
D = zeros(length(X),k);

for i = 1:k
    % Ellipsoid radius
    R = (det(Me(:,:,i)))^(1/3);
    % Ellipsoid covariance matrix based on eq (x-xc)'*A(x-xc)
    A = (inv(Me(:,:,i)))^2;
    % Substract center offset
    Xd = X-oe(i,:);
    % Distances (Normalized Mahalanobis) - vectorized to work quickly,
    % should use A' but A is symmetric anyways
    D(:,i) = R*sqrt(sum(Xd.*(Xd*A),2));
end

% Find closest cluster based on M distance to center
[~,idxClust] = min(D,[],2);
% Assign labels to voxels in label image
imgLab(sub2ind(sizeIm,X(:,1),X(:,2),X(:,3))) = idxClust;