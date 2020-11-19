function [apprxCleanPsd,approxNoisePsd,alphaApprx,stopPar]=finding_min_als(Sreal,Eps,maxIter)
% ALS method for RPSD factorization
% 
% Approximate Clean and Noise PSD and the particle location vector alpha.
% 
%
% Amitay Eldar, Dec 2017
% 
% Input parameters:
% Sreal        PSD matrix to be factorized. 
% Eps          convergence term.
% maxIter      maximum iterations allowed.
% 
% Output parameters:
% apprxCleanPsd  Approximated clean PSD.
% apprxNoisePsd  Approximated noise PSD.
% alphaApprx     Particle location vector alpha.
% stopPar        stop algorithm if an error occured 

% initial parameters
sz=size(Sreal);
patchNum=sz(2);
One=ones(1,patchNum);

% % Initial vectors
SrealNormInf = vecnorm(Sreal,Inf);
[~,maxCol] = max(SrealNormInf); [~,minCol]=min(SrealNormInf);
cleanSigTmp = abs(Sreal(:,maxCol)-Sreal(:,minCol));
SrealNormOne = vecnorm(Sreal,1);
[~,minCol]=min(SrealNormOne);
noiseSigTmp = abs(Sreal(:,minCol));
S = Sreal-noiseSigTmp*One;
alphaTmp = (cleanSigTmp'*S).*(1/sum(cleanSigTmp.^2));
alphaTmp = alphaTmp-((alphaTmp<0).*alphaTmp)-(alphaTmp>1).*((alphaTmp)-ones(1,patchNum));
stopPar = 0;

%% ALS iteration proccess as written in section 5
iter=1;
while stopPar == 0
    if norm(alphaTmp,1)==0 % if alpha is zero vecotr then we start again the process.
        alphaTmp = rand(size(alphaTmp));
    end
    apprxCleanPsd = (S*alphaTmp').*(1/sum(alphaTmp.^2));
    apprxCleanPsd = apprxCleanPsd-(apprxCleanPsd<0).*apprxCleanPsd;
    S=Sreal-apprxCleanPsd*alphaTmp;
    approxNoisePsd = (S*ones(patchNum,1)).*(1/patchNum);
    approxNoisePsd = approxNoisePsd-(approxNoisePsd<0).*approxNoisePsd;
    S = Sreal-approxNoisePsd*One;
    if norm(apprxCleanPsd,1)==0 % if apprxCleanPsd is zero vecotr then we start again the process.
        apprxCleanPsd = rand(size(apprxCleanPsd));
    end
    alphaApprx = (apprxCleanPsd'*S).*(1/sum(apprxCleanPsd.^2));
    alphaApprx = alphaApprx-((alphaApprx<0).*alphaApprx)-(alphaApprx>1).*((alphaApprx)-ones(1,patchNum));
    if norm(noiseSigTmp-approxNoisePsd)/norm(approxNoisePsd)<Eps
        if norm(cleanSigTmp-apprxCleanPsd)/norm(apprxCleanPsd)<Eps
            if norm(alphaApprx-alphaTmp)/norm(alphaApprx)<Eps
                break
            end
        end
    end
    noiseSigTmp = approxNoisePsd;
    alphaTmp = alphaApprx;
    cleanSigTmp = apprxCleanPsd;
    iter = iter+1;
    if iter>maxIter
        stopPar=1;
        break
    end 
end


