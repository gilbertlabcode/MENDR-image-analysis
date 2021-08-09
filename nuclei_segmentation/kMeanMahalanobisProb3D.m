function [Dmean,Me,oe,idx] = kMeanMahalanobisProb3D(a, k, itmax, method, idx)

% kMeanMahalanobisProb3D clusters a point set into k ellipsoids based on a
% modification of the kmeans algorithm presented in "Multiple ellipse
% fitting by center-based clustering" by Marosevic and Scitovski 2015
% (https://hrcak.srce.hr/ojs/index.php/crorr/article/view/2776). This
% version of the code uses likelihood-like functions to exclude outlier
% points when fitting ellipsoids for re-clustering.
%
% USAGE:
%
%   [Dmin,Me,oe,idx,P] = kMeanMahalanobisProb3D(a, k, itmax, method, idx)
%
% INPUTS:
%
%   a:          Matrix (nx3) with coordinates of point cloud to be clustered
%   k:          Number of clusters (ellipsoids) to fit
%   itmax:      Maximum number of iterations
%   method:     Method for fitting ellipsoid. Must be 'SOD' or 'HES'
%   idx:        Initial cluster index of points in a (starting values)
%
% OUTPUTS:
%
%   Dmean:      Mean distance of all points to their ellipsoid centers
%   Me:         3x3xk matrix that contains ellipsoid shape parameters for
%               each of the k ellipsoids
%   oe:         3xk matrix that contains ellipsoid centroids for each k
%               ellipsoids
%   idx:        Index vector indicating to which cluster each point belongs
%               to
%
% AUTHOR:
%   Jose L. Cadavid, University of Toronto, 2021

% Initialize cell arrays for the parameters of the ellipsoids
S = cell(k,1);
[np,~] = size(a);
% Initialize algebraic distances for clustering (matrix with np x k
% elements)
D = zeros(np,k);
% Initialize likelihoods
L = zeros(np,k);
P = zeros(np,k);

% Loop for itmax iterations. error is a variable to track how the
% clustering assignment changes (to check for convergence)
it = 1;
error = 100;


while (it <= itmax) && (error > 0)
    
    %% Step 1: Update step - get parameters of ellipse centers for the given clusters
    
    % Update parameters for k ellipsoids: Use direct fit method and get
    % structure with parameters for each ellipsoid
    for i = 1:k
        % Only refit if new cluster has more than 40 points to prevent
        % poorly conditioned fits that might be transient. This number is
        % pretty arbitrary for now and might be tuned or removed if
        % exceptions are handled differently. It's just there because an
        % ellipsoid needs 9 points at least to get parameters
        if sum(idx==i)>40
            S{i}=hyperellipsoidfit(a(idx==i,:),'auto',method);
        end
    end
    
    %% Step 2: Assignment - assign points to closest ellipse center
    
    % Get algebraic distance from all points to each ellipsoid center. Using
    % the regressors (Design matrix) and the parameters from the ellipsoid
    % fit (not normalizing data)
    R = zeros(k,1);
    for i = 1:k
        % Ellipsoid radius
        R(i) = (det(S{i}.Me))^(1/3);
        % Ellipsoid covariance matrix based on eq (x-xc)'*A(x-xc) and the
        % sphere to ellipsoid map Me
        Sig = (inv(S{i}.Me))^2;
        % Substract center offset from points of ellipsoid
        Xd = a - S{i}.oe';
        % Get squared algebraic distances normalized by the squared radius
        % of ellipsoid
        D(:,i) = (R(i)^2-R(i)^2*(sum(Xd.*(Xd*Sig),2))).^2;
        
        % Get probabilities: assuming algebraic distances follow exponential
        % distribution (they seem to!)
        
        % Rate parameter estimated for points in each cluster
        mu = 1./mean(D(idx==i,i));
        % Likelihoods
        L(:,i) = mu*exp(-mu*D(:,i));
    end
    
    % Probabilities that a voxel belongs to a cluster: unweighted ratio of
    % likelihoods
    for i = 1:k
        P(:,i) = L(:,i)./sum(L,2); 
    end

    % Assign points to ellipsoid only if the probability of that point being
    % part of that ellipsoid is high (normalized likelihood). This will
    % typically happen for the ellipsoid center closest to the point.
    % However, in some situations the point is close to more than one
    % ellipsoid (i.e. undecided point) so it is not clustered to prevent it
    % from affecting the fit of the ellipsoids in the next round.
    
    idx2=k*ones(np,1)+1;% cluster k+1 is a void cluster
    Plim = 0.99;
    for i = 1:k
        %Assign based on probability
        idx2(P(:,i)>=Plim)=i;
    end
    
    % Count the number of points that are assigned to a different cluster.
    % We stop when no assignment changes
    error = sum(idx~=idx2);
    % Update
    idx = idx2;
    it = it+1;
end

% Add the closest distances by number of points
[Dmin,~] = min(D,[],2);
Dmean = sum(Dmin)/length(a);

% Extract the mapping matrices and offsets from the best ellipsoids. Store
% the mapping matrices as a 3*3*k matrix, and oe as a k*3 matrix
Me = zeros(3,3,k);
oe = zeros(k,3);
for i = 1:k
    Me(:,:,i)=S{i}.Me;
    oe(i,:)=S{i}.oe';
end

