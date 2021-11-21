function [phi,mu0_est,cov0_est,mu1_est,cov1_est] = ASOCEM(I0,downscale_size,area_mat_sz,contamination_criterion,fast_flag,maxIter)
% This function segment contamination in cryo_EM micrographs

% INPUT
% I0               micrograph image .mrc format
% area_mat_sz      area for covarience estimation   
% maxIter          maximum number of iteration allowed

% OUTPUT
% phi          chan vese level set function. phi>0 is area zero and phi<=0 is area 1


%% intialize segmentation parameters
dt = 10^(0); % time step
mu = 10^(-2); %curve length regularization
nu = 0; %curve area regularization
eta = 10^(-8); %prevent zero curvature (division by 0)
Eps = 1; %epsilon (for delta function)
tol = 10^(-3);
% maxImgSz = 200; % downsampling parameter
if mod(area_mat_sz,2)==0 %% cov_mat_sz have to be odd so cov matrix size will be odd
    area_mat_sz = area_mat_sz-1;
end
cov_mat_sz = area_mat_sz^2;

%% preprocess image
% size rescale
szScaling = downscale_size/max(size(I0));
I = cryo_downsample(I0,floor(szScaling*size(I0)));
% image size should be odd
if mod(size(I,1),2)==0
    I = I(1:end-1,:);
end
if mod(size(I,2),2)==0
    I = I(:,1:end-1);
end


%% intialize phi as liphsiczh circ
[X,Y] = meshgrid(-floor(size(I,2)/2):1:floor(size(I,2)/2),-floor(size(I,1)/2):1:floor(size(I,1)/2));
phi_0 = min(size(I)/3)^2 - (X.^2 + Y.^2);
phi_0 = phi_0./max(abs(phi_0(:)));


%% chan vese time process
[phi,mu0_est,cov0_est,mu1_est,cov1_est] = chan_vese_process(I,phi_0,cov_mat_sz,dt,mu,nu,eta,Eps,maxIter,tol,fast_flag);
if phi == ones(size(phi)) % we dont want to use this micrograph
    return
end

if contamination_criterion == 0 % the smaller area will considered to be the contamination
    size_0 = sum(phi>0,'all');
    size_1 = sum(phi<0,'all');
    if size_0 >= size_1 % changing areas 
        phi = -1*phi;
    end
else % the lower mean will be considerd as contamination
    tmp_0 = I(phi>0);
    m_0 = mean(tmp_0(:));
    tmp_1 = I(phi>0);
    m_1 = mean(tmp_1(:));
    if m_0 >= m_1 % changing areas 
        phi = -1*phi;
    end  
end
end

