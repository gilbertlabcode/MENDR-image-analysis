function [param, tOpt, flag] = fitCDF(centers, counts, plt)

% fitCDF takes histogram counts of an image and attempts to reproduce the histogram
% as a mixture of two log-normal distributions and returns their parameters and an
% optimual threshold for separating both classes. We also return a flag indicating
% if the segmentation is reliable.
%
% USAGE:
%   [param, tOpt, flag] = fitCDF(centers, counts, plt)
%
% INPUTS:
%
%   centers:  Vector with bin centers of the histogram of image to be thresholded
%   counts:   Vector of counts of bins corresponding to "centers"
%   plt:      Optional argument (boolean) indicating if we should make plots
%
% OUTPUTS:
%
%   param:    Vector of parameters of the fitted log-normal distributions
%             [mu_f,sigma_f, mu_b, sigma_b, theta]
%   tOpt:     Optimal threshold estimated based on maximum class uncertainty
%   flag:     Boolean that indicates if the thresholding is not reliable (1)
%
% NOTES:
%
% This method does nto use the more rigorous expectation-maximization method but
% we found it works quite well with 4x magnification images of our engineered muscle
% tissues.
%
% AUTHOR:
%   Jose L. Cadavid, University of Toronto, 2021

%% Solve optimization

% Order of parameters [mu_f,sigma_f, mu_b, sigma_b, theta]
% mu and sigma are parameters of the log-normal distribution (f and b are
% foreground and background) and theta is the mixing of both classes.

% lower bounds
lb = [0.01, 0.001, 0, 0.001, 0]';
% Upper bounds: Assuming maximum value of 255 for 8-bit images
ub = [log(255), 10, log(30), 10, 1]';
% Solve problem
options = optimoptions('fmincon','Algorithm','sqp','Display','off');

% Initialize parameters: Initial guess.
% This problem seems to be very sensitive on the initial value chosen for
% theta, with some convergence issues. To overcome this, we solve if for
% a range of initial values of theta and take the solutions that lead to the
% optimal error function

theta0 = linspace(0,1,11);
paramMat = zeros(5,11);
fvalMat = zeros(11,1);

for i = 1:11
    p0 = [log(50), 0.5, log(5), 0.5, theta0(i)]';
    [param,fval] = fmincon(@(p)errorCDF(p,centers,counts),p0,[],[],[],[],lb,ub,[],options);
    paramMat(:,i) = param;
    fvalMat(i) = fval;
end

%Look for optimal solution and use it
[~, iOpt] = min(fvalMat);
param = paramMat(:,iOpt);

%% Get optimal threshold

% Get cut point where there is maximum class uncertainty (where conditional
% probabilities are 0.5). Do by brute force starting from pixel value 255
% and going down until theta*pf(x) <0.5, sure to converge to the top
% solution

% Parse optimal parameters
mu_f=param(1);
sigma_f=param(2);
mu_b=param(3);
sigma_b=param(4);
theta=param(5);

% Calculate class uncertainty based on the partial distributions of foreground
% and background as modelled with the fitted parameters
pf = @(g) lognpdf(g,mu_f,sigma_f);
pb = @(g) lognpdf(g,mu_b,sigma_b);
p = @(g) theta*pf(g)+(1-theta)*pb(g);
h = @(g) -theta*pf(g)./p(g).*log(theta*pf(g)./p(g))-(1-theta)*pb(g)./p(g).*log((1-theta)*pb(g)./p(g));

% Loop through thresholds to find peak
tOpt = 255;
pOpt = theta*pf(tOpt)/p(tOpt);
while pOpt>0.5
    tOpt = tOpt-1;
    pOpt = theta*pf(tOpt)/p(tOpt);
end

% Check whether the optimal threshold is between the medians
% of the two populations. We initially were verifying the monotonicity
% of the partial distributions but this criterion was too stringent
m1 = exp(param(1));
m2 = exp(param(3));
% If the optimal threshold is NOT between the two medians, we raise a flag to
% indicate the classes are not well separated (i.e. the threshold is not reliable)
flag = (tOpt<=min(m1,m2)) || (tOpt>=max(m1,m2));

%% Plot: Optional if plt arg=1

% Set plot to 0 by default if argument is not included
if nargin == 1
    plt = 0;
end

if plt
    subplot(1,2,1)
    hold on
    plot(centers,theta*lognpdf(centers,mu_f,sigma_f),centers,(1-theta)*lognpdf(centers,mu_b,sigma_b))

    subplot(1,2,2)
    plot(centers,theta*lognpdf(centers,mu_f,sigma_f)./(theta*lognpdf(centers,mu_f,sigma_f)+(1-theta)*lognpdf(centers,mu_b,sigma_b)),'r','linewidth',1.2);
    hold on
    plot(centers,(1-theta)*lognpdf(centers,mu_b,sigma_b)./(theta*lognpdf(centers,mu_f,sigma_f)+(1-theta)*lognpdf(centers,mu_b,sigma_b)),'b','linewidth',1.2);
    plot(centers,h(centers),'k','linewidth',1.2)
    legend('Foreground','Background','Class uncertainty');
end

%% Function to minimize: Least squares formulations to the CDF. Less elegant than maximum likelihood but easier to formulate.
function e = errorCDF(p,x,y)

%Get parameters
mu_f = p(1);
sigma_f = p(2);
mu_b = p(3);
sigma_b = p(4);
theta = p(5);
%Calculate estimated CDF
CDF = theta*logncdf(x,mu_f,sigma_f)+(1-theta)*logncdf(x,mu_b,sigma_b);
%Get error sum (least squares)
e=dot((CDF-y),(CDF-y));
