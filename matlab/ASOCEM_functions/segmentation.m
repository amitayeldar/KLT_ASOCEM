function [I,phi_0,phi,mu0_est,R0_est,mu1_est,R1_est] = CVS_EM(I0,cov_mat_sz,maxImgSz,smoothing,maxIter,tol)

%% intialize segmentation parameters
dt = 0.1; % time step
mu = 1; %curve length regularization
nu = 1; %curve area regularization
eta = 10^(-8); %prevent zero curvature (division by 0)
Eps = 1; %epsilon (for delta function)

%% preprocess image
% size rescale
maxSz = max(size(I0));
szScaling = maxImgSz/maxSz;
I = cryo_downsample(I0,floor(szScaling*size(I0)));
% image size should be odd
if mod(size(I,1),2)==0
    I = I(1:end-1,:);
end
if mod(size(I,2),2)==0
    I = I(:,1:end-1);
end

% smoothing
if smoothing==1
%     I = imgaussfilt(I,1);
    I = imfilter(I,ones(3,3)/3^2);
%     I = mat2gray(I);
%     I = imadjust(I);
%      I = histeq(I);
%     I = adapthisteq(I);
end

%% intialize phi as liphsiczh circ
[X,Y] = meshgrid(-floor(size(I,2)/2):1:floor(size(I,2)/2),-floor(size(I,1)/2):1:floor(size(I,1)/2));
phi_0 = min(size(I)/3)^2 - (X.^2 + Y.^2);
phi_0 = phi_0./norm(phi_0,'fro');


%% chan vese script
[phi,mu0_est,R0_est,mu1_est,R1_est] = chan_vese_process(I,phi_0,cov_mat_sz,dt,mu,nu,eta,Eps,maxIter,tol);


end

