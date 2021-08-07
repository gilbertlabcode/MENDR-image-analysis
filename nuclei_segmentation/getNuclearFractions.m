function NucStats=getNuclearFractions(Nucn,GFPn,SAAn,scale)

%Nuc, GFP and SAA are the name of the tif stacks to be read. Must have been
%resliced in imageJ to become isotropic

%%% Read images and convert to nuclei to double but keep GFP and SAA in
%%% 8-bit
Nuc=double(readTiffStack(Nucn));
GFP=uint8(readTiffStack(GFPn));
SAA=uint8(readTiffStack(SAAn));

%%% Binarize and segment nuclei with the iterative k-cluster
[labelIm,dNuc,oeNuc]=segmentNucleiEllipsoid(Nuc,scale);
%Get stats of nuclei (only volume)
NucStats=regionprops3(labelIm,'Volume','VoxelIdxList');
%%% Remove any nuclei with zero volume that might have been mislabeled
NucStats(NucStats.Volume(:)==0,:)=[];

%%% Heavily filter GFP and SAA with a median filter
GFPf=medfilt3(GFP,[5 5 5]);
SAAf=medfilt3(SAA,[5 5 5]);

%%% Get an adequate threshold for GFP ans SAA (for now using a matlab
%%% implementation o the triangle method which is less conservative than
%%% Otsu's method)

%Get histograms to feed the triangle method
[GFPhist,~]=imhist(GFPf);
[SAAhist,~]=imhist(SAAf);

%Get threshold
tGFP=triangle_th(GFPhist,256)*max(GFPf(:));
tSAA=triangle_th(SAAhist,256)*max(SAAf(:));

%%% To define whether a nucleus is inside the fiber we say that at least
%%% 75% of its volume must be masked in the GFP or SAA channels (i.e. 75%
%%% of the pixels in the nucleus must be brighter than the thresholds).
%%% This is to allow for small holes that appear during threshold or
%%% imperfect color alignment from the confocal and can be tuned to be more
%%% or less stringent

for i=1:height(NucStats)
    region=NucStats.VoxelIdxList{i,1};
    Vol=NucStats.Volume(i);
    %Number of pixels in the region brighter than the threshold
    nGFP=sum(GFPf(region)>tGFP);
    %Number of pixels in the region brighter than the threshold
    nSAA=sum(SAAf(region)>tSAA);
    
    %Is the nucleus considered GFP or SAA?
    NucStats.GFP(i)=(nGFP/Vol>=0.75);
    NucStats.SAA(i)=(nSAA/Vol>=0.75);
    NucStats.SAAGFP(i)=(nGFP/Vol>=0.75)&&(nSAA/Vol>=0.75);
end


% %%% Get the mean SAA and GFP intensity in the region comprised by each
% %%% nucleus
% for i=1:height(NucStats)
% mGFP(i)=mean(GFPf(NucStats.VoxelIdxList{i,1}));
% mSAA(i)=mean(SAAf(NucStats.VoxelIdxList{i,1}));
% end
% 
% %%% Define cuts for GFP and SAA being positive by getting the Otsu's
% %%% threshold for the mean value of the nuclei: This can be tuned by
% %%% selecting a lower multiple of the threshold if the method is
% %%% underestimating the fractions
% 
% SAAthr=graythresh(mSAA/max(mSAA(:)))*max(mSAA(:));
% GFPthr=graythresh(mGFP/max(mGFP(:)))*max(mGFP(:));
% 
% %%% Total counts
% res.nTotal=height(NucStats);
% res.SAApos=sum(mSAA(:)>SAAthr);
% res.GFPpos=sum(mSAA(:)>GFPthr);
% res.SAAGFPpos=sum((mSAA(:)>SAAthr).*(mGFP(:)>GFPthr));