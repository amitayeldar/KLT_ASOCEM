function label = clusterVariance(I,k)
%%
% function creates an intial guess according to image variance. the image
% is filtered to create an standard deviation image. the hidtogram of the
% std image is found and fitted for a gaussian. using the gaussian the std
% image is thresholded to create the intial guess of I

%%
if ~exist('k','var')
    k = 5;
end

% std image
s = movingStd2(I,k);

% histogram of std image
[y,x] = histcounts(s);
x = x(2:end);

% fit histogram to gaussian to identify the background
gauss = fit(x',y','gauss1');
% define probabilty function
P = @(x) gauss.a1*exp(-((x-gauss.b1)/gauss.c1).^2);

% set std threshold to expected value of gaussian + 1 standard deviation
label = (P(s) < 0.5);

end

