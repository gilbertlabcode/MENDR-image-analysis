function idxG=getUniquekMeans(a,k,nTest)

%%% Repeat kmeans clustering nTest times on the data in matrix a for k clusters
%%% and return the indices of non-repeated clusters: kmeans might converge
%%% to different solutions and they affect the ellipsoid clustering, so we
%%% can prefilter them to make sure we only test significantly different
%%% initial conditions. This saves lots of time since sometimes there is a
%%% unique kmeans solution

%Run kmeans once: This cluster is accepted for sure and it is a reference.
%idxG is a matrix where columns represent the indices of clustering given
%by the unique kmeans cluster runs. sumD is a vector of the sum of
%distances from points in a cluster to its center. If we average those for
%each cluster it is equivalent to a mean cluster radius

[idxG(:,1),~,sumD]=kmeans(a,k,'MaxIter',200);

%Get the number of points in each cluster
N=zeros(k,1);
for j=1:k
    N(j,1)=sum(idxG(:,1)==j);
end

%Reference distance as sum of all distances divided by number of points -
%mean total radius
sumDRef=sum(sumD)/length(a(:,1));

C=mean(a,1);
sumDRef=mean(sqrt(sum((C-a).^2,2)),1);

%Divide sum of distances to centers by number of points per cluster; this
%meands that now the distances are equivalent to a mean cluster radius
sumD=sumD./N;
%Sort these distances. Sometimes kmeans returns the same solution but in
%different order, so we sort these distances to compare them consistently
sumDG(:,1)=sort(sumD);

%Globel index for storing
ii=2;

%Attempt to get other different clusters
for i=2:nTest
    [idx,~,sumDummy]=kmeans(a,k,'MaxIter',200);
    
    %Get the number of points in each cluster
    for j=1:k
        N(j,1)=sum(idx==j);
    end
    %Divide sum of distances to centers by number of points per cluster; this
    %meands that now the distances are equivalent to a mean cluster radius
    sumDummy=sumDummy./N;
    
    %Sort sum of residuals to always compare the same clusters (order can
    %change, but values shouldn't)
    sumDummy=sort(sumDummy);
    %Get difference between clusters relative to total sum of squares
    D=sqrt(sum((sumDummy'-sumDG').^2,2))/sumDRef;
    
    %If the distance is below a certain value, we accept the new cluster as a
    %different starting point that is valid (different enough to be worth
    %testing). We can be liberal here.
    
    if all(D>0.05) %if none within 5% of any other cluster, good
        idxG(:,ii)=idx;
        sumDG(:,ii)=sumDummy;
        ii=ii+1;
    end
end
