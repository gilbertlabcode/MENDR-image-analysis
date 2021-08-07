function [Dmin,Me,oe,idx,P]=kMeanMahalanobisProb3D(a,k,itmax,plt,method,idx)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Cluster point set a=(x,y) into k ellipses based on the modified k means
%%% algorithm with Mahalanobis distance from paper "Multiple ellipse
%%% fitting by center-based clustering "
%method must be either 'HES' or 'SOD'

%Initialize cell arrays for the parameters of the ellipses
S=cell(k,1);

%Initialize algebraic distances for clustering (matrix with np x k
%elements)
D=zeros(length(a),k);
%Initialize likelihoods
L=zeros(length(a),k);
P=zeros(length(a),k);
[np,~]=size(a);
%%% Loop for itmax iterations
it=1;
E=100;


while it<=itmax && E>0
    
    %% Step 1: Update step - get parameters of ellipse centers for the given clusters
    
    %Update parameters for k ellipsoids: Use direct fit method and get
    %structure with parameters for each ellipsoid
    for i=1:k
        %Try fitting ellipsoid without regularization; if it doesn't work,
        %try again with small regularization.
         
        %Only refit if new cluster has more than 40 points to prevent
        %poorly conditioned fits that might be transient. This number is
        %pretty arbitrary for now and might be tuned or removed if
        %exceptions are handled differently. It's just there because an
        %ellipsoid needs 9 points at least to get parameters
        if sum(idx==i)>40
           
        S{i}=hyperellipsoidfit(a(idx==i,:),'auto',method);
      
        end
        
       
%         if ~S{i}.success
%             'ole'
%             k
%         end
    end
  
    
    
    %% Step 2: Assignment - assign points to closest ellipse center
    
    %Get algebraic distance from all points to each ellipsoid center. Using
    %the regressors (Design matrix) and the parameters from the ellipsoid
    %fit (not normalizing data)
    
    R=zeros(k,1);
    for i=1:k
        %Ellipsoid radius
        R(i)=(det(S{i}.Me))^(1/3);
        %Ellipsoid covariance matrix based on eq (x-xc)'*A(x-xc) and the
        %sphere to ellipsoid map Me
        Sig=(inv(S{i}.Me))^2;
        %Substract center offset from points of ellipsoid
        Xd=a-S{i}.oe';
        %Get squared algebraic distances normalized by the squared radius
        %of ellipsoid
        D(:,i)=(R(i)^2-R(i)^2*(sum(Xd.*(Xd*Sig),2))).^2;
        
        %Get probabilities: assuming algebraic distances follow exponential
        %distribution (they seem to!)
        
        %Rate parameter estimated for points in each cluster
        mu=1./mean(D(idx==i,i));
        %Likelihoods
        L(:,i)=mu*exp(-mu*D(:,i));

    end
    
    %Normalize weights
%     w=w/sum(w);
    
    %Probabilities that a voxel belongs to a cluster: unweighted ratio of
    %likelihoods
    for i=1:k
        P(:,i)=L(:,i)./sum(L,2); %L*w'
    end
    

    %Assign points to a cluster by finding the ellipse center they are closest
    %to (minimum algebraic distance)
%     [Dmin,idx2] = min(D,[],2);
    
    %Assign points to ellipsoid only if the probability of that point being
    %part of that ellipsoid is big (normalized likelihood)
    idx2=k*ones(np,1)+1;% cluster k+1 is a void cluster
    Plim=0.99;%0.96-(0.96-0.5)*(it-1)/itmax;
    for i=1:k
        %Assign based on probability
        idx2(P(:,i)>=Plim)=i;
%         %Exclude points that are very far away from the centroid of the point
%         %cloud (to get rid of loner points that skew the ellipsoid fit)
%         id=find(idx2==i);
%         Dc=sqrt(sum((C(:,i)'-a(id,:)).^2,2));
%         %If this distance is larger than 1.1Rellipsoid, consider outlier
%         threshold=10;
%          idx2(id(Dc>=threshold))=k+1;


        %Exclude outliers: These are sets of points that once removed will
        %change the volume of the ellipsoid by a large number. Given that
        %the algebraic distance is a power 4 of distance, even small clouds
        %of points have a high pull. To test this, we get the volume of
        %ellipsoids that are fitted to the data (excluding points within a
        %certain distance of the point being tested).
%         id=find(idx2==i);
%         Stest=hyperellipsoidfit(a(idx2==i,:),'auto',method);
%         Vm=det(Stest.Me);
%         Vtest=zeros(length(id),1);
%         for j=1:length(id)
%         aDummy=a(id,:);
%         Dc=sqrt(sum((aDummy(j,:)-aDummy).^2,2));
%         aDummy(Dc<=10,:)=[];
%         Stest=hyperellipsoidfit(aDummy,'auto',method);
%         Vtest(j)=det(Stest.Me);
%         end
%         idx2(id(100*abs(Vtest-Vm)/Vm>=20))=k+1;
%Determine which point cloud has the most points
% id=find(idx2==i);
% [listD,nc]=getBoundaryClusters(a(id,:),1.8);
% [~,idxM]=max(nc);
% idx2(id(listD~=idxM))=k+1; %ignore points away from largest cloud
    end
    
    
    
    %Count the number of points that are assigned to a different cluster.
    %We stop when no assignment changes
    E=sum(idx~=idx2);
    %Update
    idx=idx2;
    it=it+1;
end

%Add the closest distances by number of points
[Dmin,~] = min(D,[],2);
Dmin=sum(Dmin)/length(a);

%Extract the mapping matrices and offsets from the best ellipsoids. Store
%the mapping matrices as a 3*3*k matrix, and oe as a k*3 matrix
Me=zeros(3,3,k);
oe=zeros(k,3);
for i=1:k
    Me(:,:,i)=S{i}.Me;
    oe(i,:)=S{i}.oe';
end

%% Plot

if plt
    plotEllipsoids(a,idx,Me,oe)
end

