function [apprxCleanPsd,apprxNoisePsd,noiseVarApprx,R,apprxScaling,stopPar]= rpsd_estimation(noiseMc,phi_seg,patchSz,maxIter,gpu_use)
%
% RPSD estimation
% Approximate Clean and Noise RPSD per micrograph.
% 
%
% Amitay Eldar, Dec 2017
% 
% Input parameters:
% noiseMc      micrograph input. 
% PatchSz      patch size
% maxIter      maximum iterations for findingMin function.
% 
% Output parameters:
% apprxCleanPsd  Approximated clean RPSD.
% apprxNoisePsd  Approximated noise RPSD.
% alphaApprox  Approximated noise varience.
% R              list of radii where the RPSD is sampled.
% apprxScaling   scaling factor for figs later
% stopPar        stop algorithm if an error occured 


% initial parameters
microSz = min(size(noiseMc));
m = floor(microSz/patchSz); % Number of patches on each raw and coulomn.
M = m^2;% Total number of patches for each micrographs.
L = patchSz; % Number of radial RPSD samples 
S = zeros(L,M); % matrix of RPSD samples


% for scaling the RPSD
T = 1; % Nyquist sampling rate
bandLimit = pi/T;
numOfQuad = 2^7;
[quad,nodes] = lgwt(numOfQuad,-bandLimit,bandLimit);
quad = flipud(quad); nodes = flipud(nodes);
X = repmat(quad',numOfQuad,1);
Y = repmat(quad,1,numOfQuad);
rhoMat = sqrt(X.^2 + Y.^2);
rhoMat(sqrt(X.^2 + Y.^2)>bandLimit) = 0;
[rhoSamp,~,idx] = unique(rhoMat);


%% Computing the radial RPSD of each patch
Rtmp = zeros(L,1);
mSquare = m^2;
parfor k = 1:mSquare
    row = ceil(k/m);
    col = k-(row-1)*m;
    noiseMcBlock = noiseMc(1+(row-1)*patchSz:row*patchSz,1+(col-1)*patchSz:col*patchSz);
    noiseMcBlock = noiseMcBlock - mean(noiseMcBlock(:)); % we want zero mean
    psdBlock = cryo_epsdS(noiseMcBlock,1:(patchSz)^2,floor(0.3*patchSz));
    if nnz(isnan(psdBlock))~=0 % We dont want NaN values
        'Block got NaN';
    end
    [rBlock,R] = radialavg(psdBlock,L); %Avg through Radius in the unit disk.
    % calculate var through img
    blockVar = var(noiseMcBlock(:));   
        
    % calculate var through RPSD
    psdRad = abs(triginterp(rhoSamp,R*bandLimit,rBlock));
    psdMat = reshape(psdRad(idx),numOfQuad,numOfQuad);
    varPsd = (1/(2*pi)^2)*(nodes'*psdMat*nodes);

    % scaling RPSD
    scalingPsd = blockVar/varPsd;
    rBlock = scalingPsd*rBlock;
    % setting RPSD values
    if nnz(phi_seg(1+(row-1)*patchSz:row*patchSz,1+(col-1)*patchSz:col*patchSz))==0
        S(:,k) = rBlock;
    end
    % setting sampeling values
    if k==1
        Rtmp(:,k) = R;
    end
end

%% get rid of zero columns in S
S(:,all(S == 0))=[]; % removes column if the entire column is zero

%% Finding Min argument- ALS method
R = Rtmp(:,1);
Eps = 10^(-2); % convergence term
[apprxCleanPsd,apprxNoisePsd,~,stopPar] = finding_min_als(S,Eps,maxIter);

%% Apprx noisePsd
% first apprx noiseVar
if gpu_use ==1
    stedMatGpu = stdfilt(gpuArray(noiseMc),ones(patchSz));
    stedMat = gather(stedMatGpu);
else
    stedMat = stdfilt(noiseMc,ones(patchSz));
end
varMat = stedMat.^2;
varMat(phi_seg>0) = 0;
cut = (patchSz-1)/2 +1; % we take in considiration only whole blocks
varMat = varMat(cut:end-cut,cut:end-cut);
varVec = sort(varMat(:));
varVec = varVec(varVec>0);
j = floor(0.25*length(varVec));
noiseVarApprx = mean(varVec(1:j));

% approximate varience via RPSD 
T = 1; % Nyquist sampling rate
bandLimit = pi/T;
numOfQuad = 2^7;
[quad,nodes] = lgwt(numOfQuad,-bandLimit,bandLimit);
quad = flipud(quad); nodes = flipud(nodes);
X = repmat(quad',numOfQuad,1);
Y = repmat(quad,1,numOfQuad);
rhoMat = sqrt(X.^2 + Y.^2);
rhoMat(sqrt(X.^2 + Y.^2)>bandLimit) = 0;
[rhoSamp,~,idx] = unique(rhoMat);
cleanPsdNodes = abs(triginterp(rhoSamp,R*bandLimit,apprxCleanPsd));
noisePsdNodes = abs(triginterp(rhoSamp,R*bandLimit,apprxNoisePsd));
cleanPsdMat = reshape(cleanPsdNodes(idx),numOfQuad,numOfQuad);
noisePsdMat = reshape(noisePsdNodes(idx),numOfQuad,numOfQuad);
scalingPsdApprx = ((nodes'*noisePsdMat*nodes)-(4*pi()^2)*noiseVarApprx)/(nodes'*cleanPsdMat*nodes); % using apprx var
noisePsdApprxSigma = apprxNoisePsd - scalingPsdApprx*apprxCleanPsd;% using apprx var
noisePsdApprxSigma(noisePsdApprxSigma < 0) = 0;

sMean = mean(S,2);
sMeanPsdNodes = abs(triginterp(rhoSamp,R*bandLimit,sMean));
sMeanPsdMat = reshape(sMeanPsdNodes(idx),numOfQuad,numOfQuad);
sMeanVarPsd = (1/(2*pi)^2)*(nodes'*sMeanPsdMat*nodes);
cleanVarPsd = (1/(2*pi)^2)*(nodes'*cleanPsdMat*nodes);
cleanVar = sMeanVarPsd  - noiseVarApprx;
apprxScaling = cleanVar/cleanVarPsd;

% scaling the clean and noise RPSD
apprxCleanPsd = apprxScaling*apprxCleanPsd;
apprxNoisePsd = noisePsdApprxSigma;

end








