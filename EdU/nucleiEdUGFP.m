
function [gfpT,nTot,EdUpos,nTotGFP,EdUposGFP]=nucleiEdUGFP(Nuc,GFP,EdU,scale)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Simplified workflow for EdU+ nuclei in GFP cells only. Since these
%%% cells are sparse we do not need a super accurate nuclear segmentation

%%% Jose Cadavid, March 19 2020

%%% Modified on March 25 to also count GFP blobs (cells) and EdU+ GFP cells
%%% (not nuclei-specific). Modified the second threshold for GFP to only
%%% use the threshold from the histogram if it is larger than the minimum
%%% of 20 to prevent a low threshold from being applied if there is a
%%% problem with the histogram fitting

%% Get segmentation parameters: Geometry constraints for nuclear blobs
p=segmentationParametersNuclei();

%% Filter and binarize the nuclei channel

%%% Blur: Median filter for speed, no need for anisotropic diffusion since
%%% nuclei of interest are sparse
Nuc=double(medfilt3(Nuc,[3 3 3]));
%%% Binarize
NucBin=imbinarize(Nuc/255);

%%% Label connected components and remove small blobs (of less than Vmin)
%%% Tuning point, volume of blobs to be removed can be changed

%Get volume of blobs
stats = regionprops3(bwlabeln(NucBin),'Volume','VoxelIdxList');
%Transform that volume from voxels to cubic microns with the scale provided
stats.Volume=stats.Volume*scale^3;

for i=1:height(stats)
    if stats.Volume(i)<p.vMin
        %Small component, delete voxels
        NucBin(stats.VoxelIdxList{i,1})=0;
    end
end

%At this point we have a clean nuclear mask. We do not attempt to identify
%single nuclei yet.

%% Filter and binarize the GFP channel

GFP=double(medfilt3(GFP,[3,3,3]));

%Since the GFP fibers are sparse, thresholding methods might identify a lot
%of background. To get around this, I propose using Otsu's method to get an
%initial idea where the fibres are (or any fixed threshold, 20 seems to be
%a good value).

%Then, the binary mask is dilated 2 voxels which captures a lot of
%surrounding voxels (approximately doubles the volume of the mask). Then,
%we can use the histogram decomposition technique used for the 4x images to
%find a new threshold to rebinarize

%Otsu's - some images are almost empty so we make sure that the threshold
%is at least 20 (determined based on images)
gfpT=max(20,255*graythresh(GFP/255));
GBin=GFP>gfpT;

if gfpT>20 %If the threshold was set to 15 it is likely that there are no fibers here, so do not rebinarize to avoid error

    %Update March 26 2020: Deleted this section that used histogram fitting
    %to expand the GFP mask since it was unreliable. Switched to another
    %approach using voxel intensity and distance to blobs to allow
    %inclusion of voxels originally considered black. It also doesn't create a
    %HARD threshold
    
    %Dummy image from dilation
%     imDummy=imdilate(GBin,strel('sphere',2));
%     
%     %Get new threshold
%     [~, tOpt,~,~]=thresholdImage(GFP(imDummy(:)==1),0);
%     
%     %Rebinarize GFP: As long as the new threshold is still larger than 20
%     if tOpt>20
%         GBin=GFP>tOpt;
%     end
%     
%Get distance map: Distances to all black voxels to ANY surface of the GFP
%blobs (much faster than getting the distance to each blob and might be
%good enough here)
imD=bwdist(GBin);

%Get penalties: Based on normal distribution, square difference of GFP
%intensity vs mean GFP intensity of initial mask, normalized by variance of GFP
%in mask. This might cut some of the brightest pixels but we will include
%them again. There is another penalty based on the distance map and an
%normal distribution normalized by 5 pixel distance (tunable)

P=exp(-(GFP-mean(GFP(GFP(:)>gfpT))).^2/(2*var(GFP(GFP(:)>gfpT)))).*exp(-imD.^2/5^2);

%Apply threshold based on penalty (the larger the more stringent), 0.2
%seems fine.
GBin=P>0.2 | GFP>gfpT;
    
    
    
end

%Remove GFP things smaller than about a nucleus
stats=regionprops3(bwlabeln(GBin),'Volume','VoxelIdxList');
stats.Volume=stats.Volume*scale^3;

%Sanity check, sometimes there is ONE voxel that passes and then it messes up the bracket indexing for the voxel list
if all(stats.Volume<=p.vMin)
    GBin=zeros(size(GBin));
else
    for i=1:height(stats)
        GBin(stats.VoxelIdxList{i})=stats.Volume(i)>p.vMin;
    end
end

%% Mask the binarized nuclei with the GFP
NGFP=GBin.*NucBin;

%Remove any spurious chunk of nuclei that might have been included by
%eliminating tiny chunks

stats=regionprops3(bwlabeln(NGFP),'Volume','VoxelIdxList');
stats.Volume=stats.Volume*scale^3;

%Sanity check, sometimes there is ONE voxel that passes and then it messes up the bracket indexing for the voxel list
if all(stats.Volume<=p.vMin) 
    NGFP=zeros(size(NGFP));
else
    
    for i=1:height(stats)
        NGFP(stats.VoxelIdxList{i})=stats.Volume(i)>p.vMin;
    end
    
end
%At this stage overlaps in nuclei are unlikely since they are sparse, so we
%don't split blobs or anything

%% Determine EdU+ nuclei

%Filter EdU
EdU=double(medfilt3(EdU,[3 3 3]));

%To determine the nuclei that are EdU+, we get the mean EdU value in the
%region defined by each GFP+ nucleus. I've noticed paper fibers are visible
%in the EdU channel at a level between background and actual EdU.

%To determine a value to define a nucleus as EdU+ we could use a fix value
%(say 10) but we can define the background level of EdU by analysing its
%intensity in the region where no nuclei are detected! We can build a
%histogram of EdU pixels in nuclei- regions and define a threshold as the
%value where the cumulative distribution of those pixels is 97.5% or so

h=EdU(~NucBin(:));
N=sum(~NucBin,'all');

tEdU=0.5;
stop=0;
%Get that 957.5% threshold (modified on March 31, 2020. Previously it was
%95%)
while ~stop
    stop=sum(h<tEdU)/N>0.975;
    tEdU=tEdU+0.5;
end


% h=EdU(NucBin(:));
% tEdU=max(h)*graythresh(h/max(h))
% 
% 
% stats=regionprops3(bwlabeln(NucBin),'Volume','VoxelIdxList');
% h=[];
% for i=1:height(stats)
%     h(i)=mean(EdU(stats.VoxelIdxList{i}));
% end
% 
% tEdU=max(h)*graythresh(h/max(h))




%Get mean EdU value in GFP nuclei
stats=regionprops3(bwlabeln(NGFP),'Volume','VoxelIdxList');
EdUM=zeros(height(stats),1);

for i=1:height(stats)
    EdUM(i)=mean(EdU(stats.VoxelIdxList{i}));
end


%EdU+ percentage
EdUpos=sum(EdUM>tEdU);
nTot=height(stats);

%% Determine EdU+ GFP objects: Same as with the nuclei, but using the entire GFP blob as a mask
%Get mean EdU value in GFP blob
% stats=regionprops3(bwlabeln(GBin),'Volume','VoxelIdxList');

%Ver 5 update March 31, 2020: Getting the Mean EdU signal now on the
%gegions where there are nuclei in the blob as opposed to the mean EdU
%value in the entire blob (more precise now). Made the threshold for EdU
%tighter by increasing it to the top 2.5% intensity as opposed to the
%previous 5%

EdUM=zeros(height(stats),1);
labelGFP=bwlabeln(GBin);
nTotGFP=max(labelGFP(:));

for i=1:nTotGFP
%     EdUM(i)=mean(EdU(stats.VoxelIdxList{i}));
    EdUM(i)=mean(EdU(NGFP(:)==1 & labelGFP(:)==i))>tEdU;
end


%EdU+ percentage
% EdUposGFP=sum(EdUM>tEdU);
EdUposGFP=sum(EdUM);



% volshow(GBin);
