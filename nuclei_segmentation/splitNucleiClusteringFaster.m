function [dOpt,oeOpt,idxVox,idxClust] = splitNucleiClusteringFaster(imgBin, imgBound, statsImg, PARAM, scale)

% splitNucleiCulsteringFaster attempts to split a binary blob into k ellipsoid, 
% where different k are tested and the split for k ellipsoids is done with a modified
% Mahalanobis k-means algorithm.
%
% USAGE:
%
%   [dOpt,oeOpt,idxVox,idxClust] = splitNucleiClusteringFaster(img, imgBound, statsImg, PARAM, scale)
%
% INPUTS:
%
%   imgBin:      3D binary image of a blob to be split
%   imgBound:    3D binary image of the surface of the blob
%   statsImg:    Table with properties of the blob (volume, voxel list)
%   PARAM:       Structure with list of parameters produced by
%                segmentationParametersNuclei()
%
% OUTPUTS:
%
%   nucStats:    Table with volume of segmented nuclei. It also indicates
%                if nuclei are GFP, SAA or double-positive
%
% AUTHOR:
%   Jose L. Cadavid, University of Toronto, 2021

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% 
%%% Takes as arguments a binary image with the blob, an image of the
%%% boundary voxels of that blob and some statistics on that blob (surface
%%% area, volume, location).
%%% Jose Cadavid, 2020. University of Toronto

%% Get information from inputs

% Get the volume of the blob, the index that indicates its voxels and the
% coordinates of said voxels
Vol = statsImg.Volume;
idxVox = statsImg.VoxelIdxList{1};
[x,y,z] = ind2sub(size(imgBin), idxVox);
% Get the location of the boundary voxels. Using indices as coordinates.
% These locations are in index coordinate, not in microns. Consequently, the
% parameters of ellipsoids are also in voxel unit.
[xb, yb, zb] = ind2sub(size(imgBound), find(imgBound(:)==1));

%% Split blob into ellipsoids

% Check the ellipticity of the blob. If it is higher than 0.95 and the volume
% is smaller than the max volume, we consider the blob to be a single nucleus 
% and don't bother trying to split it. If not, we attempt to split it between 
% kmin and kmax ellipsoids (vmax and the threshold for ellipticity can be 
% tuned, but 0.95 is a very good ellipsoid)

if (Vol < PARAM.vMax) && (blobEllipticity([x,y,z],Vol,statsImg.SurfaceArea) > 0.95)
    % Blob is very much an ellipsoid and its size is consistent with the max
    % volume, so we fit a single ellipsoid directly
    S = hyperellipsoidfit([xb yb zb],'auto','SOD');
    % Get center and axes of ellipsoid based on the mapping matrix Me
    oeOpt = S.oe';
    dOpt = eig(S.Me)';
    % Create cluster ID. Everything belongs to cluster 1 since there is only
    % one ellipsoid
    idxClust=ones(length(idxVox),1);
    
else
    %% Break blob into ellipsoids

    % The blob cannot be considered a single ellipsoid and we attempt to
    % split it. The range of k to be explored is defined based on the volume
    % of the blob and the median volume of the nuclei (can be tuned depending on the expected nuclei size).
    
    % Create a window around the number of ellipsoids expected based on
    % blob size
    kmin = max(1,floor(Vol/PARAM.vMedian)-7);
    kmax = ceil(Vol/PARAM.vMedian)+2;
    
    %%% Create mesh of points for evaluating volume overlap of ellipsoids and image matching
    
    % Get box boundaries of contour points
    xmin = min(xb);
    xmax = max(xb);
    ymin = min(yb);
    ymax = max(yb);
    zmin = min(zb);
    zmax = max(zb);
    % Create mesh for those points: coordinates flipped to match the index
    % coordinate of the image
    [ym, xm, zm] = meshgrid(ymin:ymax,xmin:xmax,zmin:zmax);
    xm = xm(:);
    ym = ym(:);
    zm = zm(:);
    Xm = [xm, ym, zm];
    
    %%% Try to fit k ellipsoids to the data based on the modified kmeans with Mahalanobis distance
    
    % Initialize optimal number of ellipsoids (initially guessed as kmin)
    kOpt = kmin;
    
    for k = kmin:kmax
       
        %%% Cluster points for k ellipsoids

        % Run k means 20 times to get the unique clusters it produces
        % as opposed to running everything 20 times. This way, we can
        % run more tests and be sure we pass a diversity of clusters to
        % the algorithm, and also avoig repeating clusters (more
        % efficient). This is done since sometimes kmeans returns different
        % clusters
        idxG = getUniquekMeans([xb yb zb],k,20);
        % Number of distinct results from kmeans (within 5% of each other)
        [~,nTest] = size(idxG);
        
        % Initialize variables for storing data from clustering attempts
        MeDummy = cell(nTest);
        oeDummy = cell(nTest);
        DminDummy = zeros(1,nTest);
        fPenDummy = zeros(1,nTest);
        
        
        % The ellipsoid-clustering algorithm uses kmeans as an initial
        % approximation and might converge to sub-optimal points, therefore we
        % run it on all the unique cluster results (nTest) and choose the
        % best one
         for i = 1:nTest
            % Attempt clustering with each distinct cluster
            [DminDummy(i),MeDummy{i},oeDummy{i},~] = kMeanMahalanobisProb3D([xb yb zb],k,200,'SOD',idxG(:,i));
            % Evaluate solution based on nuclei overlap and volume of the
            % blobs: Get the penalized objective function that we want to
            % optimize
            fPenDummy(i) = evalClusterEllipsoid(imgBin,Xm,MeDummy{i},oeDummy{i},PARAM,scale);
        end
        
        
        % Find best fit for the given k ellipsoids among all initial clusters and store it in a global variable
        % (only store the parameters of the ellipsoids)
        [fPen(k-kmin+1),idxMax] = max(fPenDummy);
        Me{k-kmin+1} = MeDummy{idxMax};
        oe{k-kmin+1} = oeDummy{idxMax};
        
        % If we have already tested the first kmin, check whether the
        % penalized objective function increases. If yes, we update the
        % optimal kOpt and iterate again. If not, we stop (we reached a
        % local maximum)
        if k > kmin
            if fPen(k-kmin+1) > fPen(k-kmin) %better clustering
                kOpt = k;
            else
                break
            end
        end
        
    end % End of testing different number of clusters

    %% Return the optimal clustering and split the blob according to the ellipsoids found
    
    % Store the parameters of the ellipsoids. Not interested in the rotation of
    % ellipsoids at this point, so instead of passing the mapping matrix I can
    % pass the semiaxes of the ellipsoid. At this point it is possible that
    % some small nuclei have been segmented or have appeared (despite the penalty terms) so we can
    % remove them. The voxels will still be lumped with the closest
    % ellipsoid when splitting the blob. Only do if there is more than one
    % ellipsoid in blob!!
    MeOpt = Me{kOpt-kmin+1};
    oeOpt = oe{kOpt-kmin+1};
    dOpt = zeros(kOpt,3);
    for i = 1:kOpt
        dOpt(i,:) = eig(MeOpt(:,:,i))';
        Vol(i) = 4/3*pi*prod(dOpt(i,:))*scale^3;
    end
       
    if kOpt>1
        % Look if any ellipsoid is quite small and remove them
        idxRemove = Vol<PARAM.vMin/2;
        % To prevent errors, keep at least one ellipsoid (the largest one)
        if all(idxRemove)
            idxRemove(Vol == max(Vol)) = 0;
        end
        dOpt(idxRemove,:) = [];
        oeOpt(idxRemove,:) = [];
        MeOpt(:,:,idxRemove) = [];
    end
    
    % Label voxels in blob according to optimal ellipsoids: Return linear index
    % of voxels and their cluster label
    [~,idxClust] = labelSplitNuclei(size(imgBin),[x y z],MeOpt,oeOpt);
    
end

