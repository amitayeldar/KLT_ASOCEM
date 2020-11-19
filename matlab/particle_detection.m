function [numOfPickedPar,numOfPickedNoise] = particle_detection(noiseMc,phi_seg,eigFun,eigVal,numOfFun,noiseVar,mcSz,mgScale,radMat,mgBigSz,patchSzPickBox,patchSzFun,num_of_particles,num_of_noise_images,coordinatsPathParticle,coordinatsPathNoise,microName,thresh,gpu_use)
% Construct the scoring matrix and then use the picking_from_scoring_mat
% function to pick particles and noise images.
% 
% Amitay Eldar, Dec 2017
% 
% Input parameters:
% noiseMc               micrograph image          
% eigFun                eigenfunction as coulmns of the matrix eigFun
% eigVal                eigenvalues.
% numOfFun              number of eigenfunction to use out of eigFun.
% noiseVar              noise estimated varience
% mcSz                  size of downcampled micrograph
% mgScale               downsampling scalar
% radMat                see KLTpicker preprocess
% mgBigSz               size of micrograph
% patchSzPickBox        particle box size to use
% patchSzFun            the patch size used to construct the eigenfunctions 
% num_of_particles      amount of particle to pick. -1 to pick all.
% num_of_noise_images   amount of noise images to pick.
% boxPathParticle       path to the particle coordinate files
% boxPathNoise          path to the noise images coordinate files
% microName             name of micrograph  
% thresh                threshold for the picking out of the scoring mat.
% gpu_use               if set to 1 then use the GPU

% Output parameters:
% coordinate files in boxPathparticle and boxPathNoise.
% numOfPickedPar        number of picked particles.
% numOfPickedPar        number of picked noise images.

    eigFunStat = eigFun(:,1:numOfFun); eigValStat = eigVal(1,1:numOfFun);

    for i = 1:numOfFun
        tmpFun = reshape(eigFunStat(:,i),patchSzFun,patchSzFun);
        tmpFun(radMat>floor((patchSzFun-1)/2)) = 0; % puting zero outside the disk
        eigFunStat(:,i) = tmpFun(:);
    end

    [Q,ar] = qr(eigFunStat);
    ar = ar(1:numOfFun,1:numOfFun);
    kapa = ar*diag(eigValStat)*ar'+noiseVar*eye(numOfFun);
    kapaInv = inv(kapa);
    Tmat=(1/noiseVar)*eye(numOfFun)-kapaInv;
    [P,D] = eig(Tmat);
    Pfull = zeros(size(Q)); 
    Pfull(1:size(Tmat,1),1:1:size(Tmat,2)) = P;
    mu = logdet((1/noiseVar)*kapa);
    noiseMcGpu = gpuArray(single(noiseMc));

    % Iterating on the patches
    lastBlockRow = mcSz(1)-patchSzFun+1;
    lastBlockCol = mcSz(2)-patchSzFun+1;
    numOfPatchRow = length(1:1:lastBlockRow);
    numOfPatchCol = length(1:1:lastBlockCol);
    
    % Computes Loglikelyhood test of each patch with conv
    QP = Q*Pfull;
    if gpu_use==1
        logTestMat = gpuArray(zeros(numOfPatchRow,numOfPatchCol));
    else
        logTestMat = zeros(numOfPatchRow,numOfPatchCol);
    end

    for i = 1:numOfFun
        qptmp = reshape(QP(:,i),patchSzFun,patchSzFun); 
        qptmp = flip(flip(qptmp,1),2);
        if gpu_use==1
            scoreTmp = conv2(noiseMcGpu,qptmp,'valid');  % compute the ith component of the diagonelized bilinear form
            logTestMat = logTestMat + D(i,i)*abs(scoreTmp.^2); % D(i,i) is real.
        else
            scoreTmp = conv2(noiseMcGpu,qptmp,'valid'); % compute the ith component of the diagonelized bilinear form
            logTestMat = logTestMat + D(i,i)*abs(scoreTmp.^2);% D(i,i) is real.
        end
    end
    logTestMat = gather(logTestMat)-mu;
    % this part is for centering
    patchSzN = patchSzFun;
    if gpu_use==1
        neigh = gpuArray(ones(patchSzN));
        logTestN = gather(conv2(logTestMat,neigh,'valid'));
    else
        neigh = ones(patchSzN);
        logTestN = conv2(logTestMat,neigh,'valid'); 
    end
    [numOfPickedPar,numOfPickedNoise] = picking_from_scoring_mat(logTestN,phi_seg,patchSzFun,patchSzN,patchSzPickBox,num_of_particles,num_of_noise_images,mgBigSz,mgScale,coordinatsPathParticle,coordinatsPathNoise,microName,thresh);
end