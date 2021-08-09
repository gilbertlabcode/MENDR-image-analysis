function fPen = evalClusterEllipsoid(imgBin, Xm, Me, oe, PARAM, scale)

% evalClusterEllipsoid determines the quality of segmentation by looking at
% the fitted ellipsoids and assessing how well they fit the blob they
% describe, and penalizes this image matching based on the size of the
% ellipsoids and how much of their volume overlaps with other ellipsoids
% since real nuclei shouldn't overlap too much.
%
% USAGE:
%
%   [fPen,V,vMat,vOver,imMatch] = evalClusterEllipsoid(Im, Xm, Me, oe, p, scale)
%
% INPUTS:
%
%   imgBin:    3D binary image of blob to be split into ellipsoids
%   Xm:        Matrix describing a mesh grid
%   Me:        3x3xk matrix containing shape parameters of each k ellipsoid
%   oe:        3xk matrix containing centroids of eack of k ellipsoids
%   PARAM:     Structure containing parameters produced by
%              segmentationParametersNuclei()
%
% OUTPUTS:
%
%   fPen:      Penalized objective function (image match) to choose number
%              of ellipsoids
%
% AUTHOR:
%   Jose L. Cadavid, University of Toronto, 2021

% Initialize variables

% Number of clusters (ellipsoids)
[k,~] = size(oe);
% Volume of the individual ellipsoids
Vol = zeros(1,k);

%% Loop through ellipsoids to determine their volume and assess how they match the image

% Note: Cound be done with Monte Carlo sampling, but the mesh is good enough since
% we only need a decent approximation of volumes.

% Check which points of the proposed mesh are part of the blob in the image
inBlob = imgBin(sub2ind(size(imgBin),Xm(:,1),Xm(:,2),Xm(:,3))) == 1;
% Check if the points of the mesh are inside any ellipsoid - > (x-xc)'A(x-xc)<=1
inEllip = zeros(length(Xm),k);

for i = 1:k
    % Get the volume of each fitted ellipsoid: The fitted ellipsoid is in
    % voxels, so we multiply it by scale^3 to make the volume into microns
    Vol(i) = 4*pi*det(Me(:,:,i))/3*scale^3;
    % Ellipsoid covariance matrix based on eq (x-xc)'*A(x-xc)
    % and the sphere to ellipsoid mapping Me
    A = (inv(Me(:,:,i)))^2;
    % Substract center offset from points from mesh
    Xd = Xm - oe(i,:);
    % Tag points in the mesh that belong to the current ellipsoid
    inEllip(:,i) = double(sum(Xd.*(Xd*A),2) <= 1);
end

% With these tags, get an overlap matrix that indicates the overlap between
% pairs of ellipsoids
vMat = zeros(k*(k-1)/2,1);
ii=1;
for i = 1:k
    for j = i+1:k
        % Shared volume divided by volume of smallest blob.
        % The number of shared voxels (vShared) is multiplied by the scale^3
        % to bring it to the dimensions of actual volume
        vShared=inEllip(:,i)'*inEllip(:,j)*scale^3;
        vMat(ii)=100*vShared/min(Vol(i),Vol(j));
        ii=ii+1;
    end
end
inEllip = sum(inEllip,2);

%% Calculate constraints

% Calculate how well the ellipsoids match the image: Count the voxels that
% are in the blob AND in any ellipsoids but penalize (with pm) for voxels that are in
% tellipsoids but NOT in the blob. Normalize by the number of voxels in the
% blob (its volume). A perfect fit will be 1; this value can be negative if
% the fit is poor
imMatch = (sum(inEllip>0 & inBlob)-(PARAM.pM)*sum(xor(inEllip>0,inBlob)))/sum(inBlob>0);

% If there is only one nucleus there is no overlap (change later so I don't
% calculate this twice
if k==1
    vMat = 0;
end

%% Calculate penalties

% Penalize if the ellipsoids are outside a defined volume range. The
% parameters for the logistic equation are in the high level function
% segmentationParametersNuclei, as all other parameters

% For large ellipsoids
vo1 = PARAM.logVol(2,2);
ko1 = PARAM.logVol(2,1);
% For small ellipsoids
vo2 = PARAM.logVol(1,2);
ko2 = PARAM.logVol(1,1);
% Add all penalties: The penalty for large nuclei is not penalized as much
% as for small nuclei (which prevents oversplitting)
pV = sum(1-1./(1+exp(-ko2*(Vol-vo2))))+sum(1./(1+exp(-ko1*(Vol-vo1))));

% Penalize overlap, using logistic function: The following parameters result
% in a 5% penalty for a 10% overlap and a 95% penalty for a 20% overlap. The
% penalty is calculated for each pair of overlaps in vMat. The parameters can be tuned but
% this seems to be a reasonable choice based on the overlaps of blobs with
% two nuclei
vo = PARAM.logOver(1,2);
ko = PARAM.logOver(1,1);

pO = sum(1./(1+exp(-ko*(vMat-vo))));

% Penalize ellipsoids with high aspect ratio (to avoid very ellongated
% ellipsoids)

%Get aspect ratios
AR = zeros(k,1);
for i = 1:k
    dOpt = eig(Me(:,:,i))';
    AR(i) = max(dOpt)/min(dOpt);
end

vAR = PARAM.logAR(1,2);
kAR = PARAM.logAR(1,1);
pAR = sum(1./(1+exp(-kAR*(AR-vAR))));

% Calculate the penalized objective function: Penalize image match based on
% excess volume, or overlaps, or aspect rations
fPen = imMatch - ((PARAM.pV)*pV+(PARAM.pO)*pO+(PARAM.AR)*pAR);