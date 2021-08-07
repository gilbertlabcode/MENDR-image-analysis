function [fPen,V,vMat,vOver,imMatch]=evalClusterEllipsoid(Im,Xm,Me,oe,p,scale)

%% Initialize variables

%Number of clusters (ellipsoids)
[k,~]=size(oe);
%Volume of the individual ellipsoids
V=zeros(1,k);

%% Loop through ellipsoids to determine their volume and assess how the match the image

%Cound be done with Monte Carlo sampling, but the mesh is good enough since
%we only need a decent approximation of volumes.

%Check which points of the proposed mesh are part of the blob in the image
inBlob=Im(sub2ind(size(Im),Xm(:,1),Xm(:,2),Xm(:,3)))==1;

%Check if the points of the mesh are inside any ellipsoid - > (x-xc)'A(x-xc)<=1
inEllip=zeros(length(Xm),k);

for i=1:k
    %Get the volume of each fitted ellipsoid: The fitted ellipsoid is in
    %voxels, so we multiply it by scale^3 to make the volume into microns
    V(i)=4*pi*det(Me(:,:,i))/3*scale^3;
    %Ellipsoid covariance matrix based on eq (x-xc)'*A(x-xc)
    %and the sphere to ellipsoid mapping Me
    A=(inv(Me(:,:,i)))^2;
    %Substract center offset from points from mesh
    Xd=Xm-oe(i,:);
    %Tag points in the mesh that belong to the current ellipsoid
    inEllip(:,i)=double(sum(Xd.*(Xd*A),2)<=1);
end

%With these tags, get an overlap matrix that indicates the overlap between
%pairs of ellipsoids

vMat=zeros(k*(k-1)/2,1);

ii=1;
for i=1:k
    for j=i+1:k
        %Alternative: Shared volume divided by volume of smallest blob.
        %The number of shared voxels (vShared) is multiplied by the scale^3
        %to bring it to the dimensions of actual volume
        vShared=inEllip(:,i)'*inEllip(:,j)*scale^3;
        vMat(ii)=100*vShared/min(V(i),V(j));
        
%         vMat(ii)=100*(inEllip(:,i)'*inEllip(:,j))/(sum(inEllip(:,i)|inEllip(:,j)));
        ii=ii+1;
    end
end

inEllip=sum(inEllip,2);

%% Calculate constraints

%Calculate the % of the covered volume by the ellipoids that corresponds to
%overlaps. Voxels that are in a region of overlap have a count larger than
%1 (they are in more than one ellipsoid), whereas ALL voxels in the
%blob made of ellipsoid have a count larger than 0
vOver=sum(inEllip>1)/sum(inEllip>0)*100;

%Calculate how well the ellipsoids match the image: Count the voxels that
%are in the blob AND in any ellipsoids but penalize (with pm) for voxels that are in
%tellipsoids but NOT in the blob. Normalize by the number of voxels in the
%blob (its volume). A perfect fit will be 1; this value can be negative it
%the fit is poor
imMatch=(sum(inEllip>0 & inBlob)-(p.pM)*sum(xor(inEllip>0,inBlob)))/sum(inBlob>0);
% imMatch=(sum(inEllip>0 & inBlob)-(p.pM)*sum(inEllip>0& ~inBlob))/sum(inBlob>0); %penallize voxels that are outside of the ellipsoid
%If there is only one nucleus there is no overlap (change later so I don't
%calculate this twice
if k==1
    vOver=0;
    vMat=0;
end

%% Calculate penalties

%Penalize if the ellipsoids are outside a defined volume range. The
%parameters for the logistic equation are in the high level function
%segmentationParametersNuclei, as all other parameters

%For large ellipsoids
vo1=p.logVol(2,2);
ko1=p.logVol(2,1);
%For small ellipsoids
vo2=p.logVol(1,2);
ko2=p.logVol(1,1);
%Add all penalties: The penalty for large nuclei is not penalized as much
%as for small nuclei (which prevents oversplitting)
pV=sum(1-1./(1+exp(-ko2*(V-vo2))))+sum(1./(1+exp(-ko1*(V-vo1))));

%Penalize overlap, using logistic function: The following parameters result
%in a 5% penalty for a 10% overlap and a 95% penalty for a 20% overlap. The
%penalty is calculated for each pair of overlaps in vMat and then averaged
%as to not overpenalize partitions with more ellipsoids since the number of
%overlapping pairs increases quadratically. The parameters can be tuned but
%this seems to be a reasonable choice based on the overlaps of blobs with
%two nuclei
vo=p.logOver(1,2);
ko=p.logOver(1,1);

pO=sum(1./(1+exp(-ko*(vMat-vo))));

%Penalize ellipsoids with high aspect ratio (to avoid very ellongated
%ellipsoids)

%Get aspect ratios
AR=zeros(k,1);
for i=1:k
    dOpt=eig(Me(:,:,i))';
    AR(i)=max(dOpt)/min(dOpt);
end

vAR=p.logAR(1,2);
kAR=p.logAR(1,1);
pAR=sum(1./(1+exp(-kAR*(AR-vAR))));


%Calculate the penalized objective function: Dmin is divided by 1000 since
%that is the order of magnitude of good fit for blobs with a single
%ellipsoid (defined by ellipticity>0.95); this now makes a "good" Dmin to
%be in the order of 1. The function is penalized by nuclei overlap and
%nuclei volume, with wiegths a1 and a2

%DminV is the fit error for a number of clusters given by the volume of the
%blob divided by the median volume; this normalizes the fit of all clusters
%to what would be expected from simply taking kopt=V/Vmedian

% fPen=Dmin/DminV+a1*pV+a2*pO;
fPen=imMatch-((p.pV)*pV+(p.pO)*pO+(p.AR)*pAR);