function [param, tOpt, flag] = fitCDF(centers, counts, plt)

%Take the CDF of an image histogram and get the optimal parameters such
%that the CDF is represented as the wieghted sum of to log normal
%distributions (foreground and wbackground). The weight theta is the "true"
%(estimated) proportion of foreground in the image. In reality this
%approximation might only hold in the tails of the histogram, but it is
%useful to get priors for a more rational segmentation of the muscle images
%that we have for MENDR

%order of parameters mu_f,sigma_f, mu_b, sigma_b, theta

% thresholdImage takes in a grayscale image (8-bit) and attempts to find an optimal
% threshold to binarize it by assuming that its intensity histogram is a combination of 
% two log-normal distributions (foreground and background). The histogram is decomposed
% with a least-squares method, and the threshold is defined as the intensity at which 
% class uncertainty is maximum.
%
% USAGE:
%
%   [pArea, tOpt, flag, param] = thresholdImage(img, plt)
% 
% INPUTS:
%
%   img:    8-bit grayscale image to be thresholded
%   plt:    Optional boolean that indicates whether intermediate plots should be
%           made (1) or not (0). Default is 0
%
% OUTPUTS:
%
%   pArea:  Percentage coverage of foreground pixels
%   tOpt:   Optimal intensity threshold
%   flag:   Boolean that indicates whether the two classes of pixels are
%           not properly separated (flag = 1; threshold not reliable)
%   param:  Parameters of the fitted log-normal distributions
%
% AUTHOR:
%   Jose L. Cadavid, University of Toronto, 2021


%% Solve optimization

% order of parameters [mu_f,sigma_f, mu_b, sigma_b, theta]
%lower bounds
lb=[0.01, 0.001, 0, 0.001, 0]';
%upper bounds: Assuming maximum value of 255 for 8-bit images
ub=[log(255), 10, log(30), 10, 1]';

%Solve problem
options = optimoptions('fmincon','Algorithm','sqp','Display','off');

% Initialize parameters: Initial guess.
% This problem seems to be very sensitive on the initial value chosen for
% theta, with some convergence issues. To overcome this, we solve if for
% a range of initial values of theta and take the solutions that lead to the
% optimal error function

theta0=linspace(0,1,11);
paramMat = zeros(5,11);
fvalMat = zeros(11,1);

for i = 1:length(theta0)
    p0 = [log(50), 0.5, log(5), 0.5, theta0(i)]';
    [param,fval] = fmincon(@(p)errorCDF(p,centers,counts),p0,[],[],[],[],lb,ub,[],options);
    paramMat(:,i) = param;
    fvalMat(i) = fval;
end

%Look for optimal solution and use it
[~, iOpt] = min(fvalMat);
param = paramMat(:,iOpt);

%% Examine solution to see if histograms separate nicely

%Get cut point where there is maximum class uncertainty (where conditional
%probabilities are 0.5). Do by brute force starting from pixel value 255
%and going down until theta*pf(x) <0.5, sure to converge to the top
%solution

%Get parameters for plotting
    mu_f=param(1);
    sigma_f=param(2);
    mu_b=param(3);
    sigma_b=param(4);
    theta=param(5);
    
    %class uncertainty
    
    pf=@(g) lognpdf(g,mu_f,sigma_f);
    pb=@(g) lognpdf(g,mu_b,sigma_b);
    p=@(g) theta*pf(g)+(1-theta)*pb(g);
    h=@(g) -theta*pf(g)./p(g).*log(theta*pf(g)./p(g))-(1-theta)*pb(g)./p(g).*log((1-theta)*pb(g)./p(g));
    
tOpt=255;
p_opt=theta*pf(tOpt)/p(tOpt);
while p_opt>0.5
    tOpt=tOpt-1;
    p_opt=theta*pf(tOpt)/p(tOpt);
end

%Check whether the conditional probability of belonging to the foreground
%is monotonically increasing with pixel value. If its not, there is another
%potential point of high class uncertainty that can indicate the image does
%not have enough contrast or it has another type of object e.g. paper
%fibres that increase the apparent background

flag=0;

%Start verifying for monotonocity only from the lowest pixel value (ie.
%only check the domain of the image to avoid false flags)
i=4;%min(im(:))+1;

stop=0;
while ~stop && i<100 %&& theta*pf(i)/p(i)<0.9
   
    if (pf(i)/p(i))<(pf(i-1)/p(i-1)) %Decrease in probability with increasing pixel value, growth is not monotonial, flag and exit
        flag=1;
        stop=1;
    end
    i=i+1;
end

% x=1:255;
% ELE=theta*pf(x)./p(x);
% plot(x,ELE);
% x1=ELE(1:end-1)-0.5;
% x2=ELE(2:end)-0.5;
% S=sum(x2.*x1./(abs(x2.*x1))<=0)

% Modified flag: Check whether the optimal threshold is between the medians
% of the two populations - the original criterion was too stringent
m1=exp(param(1));
m2=exp(param(3));
flag= (tOpt<=min(m1,m2)) || (tOpt>=max(m1,m2));

%% Plot: Optional if plt arg=1

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
function e=errorCDF(p,x,y)

%Get parameters
mu_f=p(1);
sigma_f=p(2);
mu_b=p(3);
sigma_b=p(4);
theta=p(5);


%Calculate estimated CDF
CDF=theta*logncdf(x,mu_f,sigma_f)+(1-theta)*logncdf(x,mu_b,sigma_b);

%Get error sum (least squares)
e=dot((CDF-y),(CDF-y));
