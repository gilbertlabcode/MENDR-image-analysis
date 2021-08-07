function [imLabelNew,dNuc,oeNuc]=segmentNucleiEllipsoid(im,scale)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Segmentation of nuclei in 3D images based on fitting ellipsoids to
%%% binary blobs and k-means with Mahalanobis distance. Initial binarization 
%%% done by smoothing the image with a
%%% non-linear fitler to preserve edges as much as possible, and then
%%% aplying Otsu's method
%%% Written by: Jose Cadavid, 2020. University of Toronto

%%% Im is a 3D image (converted to a double matrix)

%% Get segmentation parameters
p=segmentationParametersNuclei();

%% Step 1: Binarize image by using global threshold and clean up for segmentation

%%% Blur
imF=imdiffusefilt(im,'NumberOfIterations',3,'GradientThreshold',0.2*max(im(:)));
%%% Binarize
imBin=imbinarize(imF/255);

%%% Label connected components and remove small blobs (of less than Vmin)
%%% Tuning point, volume of blobs to be removed can be changed

%Get volume of blobs
stats = regionprops3(bwlabeln(imBin),'Volume','VoxelIdxList');
%Transform that volume from voxels to cubic microns with the scale provided
stats.Volume=stats.Volume*scale^3;

for i=1:height(stats)
    if stats.Volume(i)<p.vMin
        %Small component, delete voxels
        imBin(stats.VoxelIdxList{i,1})=0;
    end
end

%%% Apply erosion and dilation to split blobs where there are small
%%% connecting regions. This also partially gets rid of tiny objects. By
%%% splitting large blobs with spherical structuring elements we can reduce
%%% computation time in clustering; especially if we split huge blobs.
%%% Small holes are also filed with morphological closing if necessary

%The size of the structuring elements will be further tweaked in a future
%iteration, and it might be calculated adaptively
imBin=splitErosion(imBin,strel('sphere',3),strel('sphere',2));

%%% Optional: Remove objects that touch border of image. For analysing the
%%% morphology of nuclei, these objects should be ignored
imBin=imclearborder(imBin>0).*imBin;

%%% At this point, the binary image is not longer binary, but the different
%%% split blobs are labelled. We go through a second round of removing tiny
%%% objects that might have appeared. We also keep a list of the blobs that
%%% are removed for relabelling

%Get volume of blobs
stats = regionprops3(imBin,'Volume','VoxelIdxList');
%Transform that volume from voxels to cubic microns with the scale provided
stats.Volume=stats.Volume*scale^3;

%Tag blobs for eliminating
list=stats.Volume<p.vMin;
imLabel=imBin; 
%Relabel without the indices of the removed blobs
 ii=1;
 for i=1:height(stats)
     if list(i) %eliminate
         imLabel(stats.VoxelIdxList{i,1})=0;
     else %relabel
         imLabel(stats.VoxelIdxList{i,1})=ii;
         %Next label
         ii=ii+1;
     end
 end
     


%%% Fill small holes of about 1/2 of the minimum nuclei size by
%%% analyzing the connectivity of the complement of the image. This ratio
%%% can be tuned or changed

% for i=1:height(stats)
%     if stats.Volume(i)<p.vMin/2 
%         %Small hole, fill voxels
%         imBin(stats.VoxelIdxList{i,1})=1;
%     end
% end

%%% Remove objects that might correspond to paper fibers. Filter by
%%% sphericity. Not removing big objects by size since there are clusters
%%% of nuclei

% stats= regionprops3(imBin,'Volume','VoxelIdxList','ConvexVolume','SurfaceArea','VoxelList');
% 
% for i=1:height(stats)
%     sph=pi^(1/3)*(6*stats.Volume(i))^(2/3)/stats.SurfaceArea(i);
%     if sph<0.7 %Observed from images, can be tuned
%         %Object with low sphericity, can be a fiber or a blob in the edge
%         %of image that looks to be cut
%         imBin(stats.VoxelIdxList{i,1})=0;
%     end
% end



%% Step 2: Split touching nuclei in the blobs with a modified k-means algorithm with Mahalanobis distance

%Relabel image: Add this point some adjacent nuclei have been labelled with
%different numbers in the image erosion algorithm, but if we relabel them
%they will be fused because they are connected (the erosion algorithm
%doesn't always create a black space between them because the dilation can
%make them touch again, which is why opening doesn't work the same way).
%Since some of the labels have been have dissappeared

%Get relevant properties of blobs (don't need convex volume or anything
%else now)
statsL=regionprops3(imLabel,'Volume','VoxelIdxList','SurfaceArea');
%Scale volume from voxels to microns
statsL.Volume=statsL.Volume*scale^3;
%Scale area from pix square to micron 2
statsL.SurfaceArea=statsL.SurfaceArea*scale^2;

%Get the boundary voxels of the blobs by convolution and counting
%neighbours. This function returns it with no labels (only 0 and 1), so we
%multiply it with the fully labelled image to label those boundary voxels
imBound=getBoundary(imLabel>0,p.nNeigh).*imLabel;

%Number of blobs
nBlobs=height(statsL);

%Initialize a new label image
imLabelNew=zeros(size(imLabel));

%Matrices for storing ellipsoid centers and semiaxes 
dNuc=[];
oeNuc=[];

%Counter for labelled nuclei
nucCount=0;

%Split blobs, one by one, based on ellipsoid clustering
for i=1:nBlobs
    
    %Get optimal split of blob in terms of ellipsoids (optimal number
    %calculated internally). Also return the center (oeOpt) and semiaxes
    %length (dOpt) of those ellipsoids. Ignoring the rotation angles but
    %could be obtained as well
    [dDummy,oeDummy,idxVox,idxClust]=splitNucleiClusteringFaster(imLabel==i,imBound==i,statsL(i,:),p,scale);
    %Get the number of new ellipsoids
    [kOpt,~]=size(oeDummy);
    %Label the voxels of the blob indicated by idxVox with the labels
    %idxClust (adding the current nuclei counter)
    imLabelNew(idxVox)=idxClust+nucCount;
    %Store semiaxes length and centers
    dNuc=[dNuc;dDummy];
    oeNuc=[oeNuc;oeDummy];
    %Update nuclei counter
    nucCount=nucCount+kOpt;
end

