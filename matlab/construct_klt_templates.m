function [eigFun,eigVal] = construct_klt_templates(rho,quadKer,quadNys,rr,sqrtrSampr,JrRho,Jsamp,cosine,sine,...
                            numOfQuadNys,maxOrder,psd,precentOfEigen,idxRsamp,gpu_use)
% Constructing the KLTpicker templates as the eigenfunctions of a given kernel.
% 
% Amitay Eldar, Dec 2017
% 
% Input parameters:
% See KLTpicker preprocess for description of the following inputs   
% rho
% quadKer
% quadNys
% rr
% sqrtrSampr
% JrRho
% Jsamp
% cosine
% sine  
% numOfQuadNys
% maxOrder         max order for the eigenfunction
% psd              particle RPSD
% precentOfEigen   how many eigenfunction to take
% gpu_use          if set to 1 then use the GPU.
% 
% Output parameters:
% eigFun           eigenfunction as coulmns of the matrix eigFun  
% eigVal           eigenvalues.


eigFunTot = zeros(numOfQuadNys,numOfQuadNys,maxOrder); % each m contain a matrix whos coulmns are eigenFun of oreder m.
eigValTot = zeros(numOfQuadNys,maxOrder); % each col contain a vector whos cordinates are the eigen values of order col.
%% finding eig fun & val
sqrt_rr = sqrt(rr);
d_rho_psd_quadKer = diag(rho).*diag(psd).*diag(quadKer);
sqrt_diag_quadNys = sqrt(diag(quadNys)); 

for N = 0:maxOrder-1
    Hnodes = sqrt_rr.*(JrRho(:,:,N+1)*(d_rho_psd_quadKer)*(JrRho(:,:,N+1))');
    tmp = sqrt_diag_quadNys*Hnodes*sqrt_diag_quadNys;
    [V,D] = eig(tmp); % find eig for similar symetric mat we will fix eig fun later  
    D = real(D);% imag val is theoriticly impossible
    [D,I] = sort(diag(D),'descend'); V = V(:,I); % sorting. note that D retruned as vector.
    D(abs(D)<eps) = 0; V(:,D==0) = 0; % this vlues can couse numeric problems;
    eigFunTot(:,:,N+1) = V; eigValTot(:,N+1) = D; % later sorting will throw away zeros.
end

%% descide how many eigs to take
% first we sort and save indexes
r_idx = [1:numOfQuadNys]'; c_idx = 1:maxOrder; 
R_idx = repmat(r_idx,1,maxOrder); C_idx = repmat(c_idx,numOfQuadNys,1);
eigValTot = eigValTot(:);   R_idx = R_idx(:);   C_idx = C_idx(:); %-check the length-%
[eigValTot,Idx] = sort(eigValTot,'descend');
R_idx = R_idx(Idx,1);   C_idx = C_idx(Idx,1);

% we sum until we get precentOfEig
sumOfEig = sum(eigValTot);
cumSumEigVal = cumsum(eigValTot/sumOfEig);
for i = 1:length(eigValTot)
    if cumSumEigVal(i)>precentOfEigen
        lastEigIdx = i;
        break
    end
end

%% Sampeling eig fun
eigVal = zeros(1,2*lastEigIdx); %if N==0 we take one fun, if N~=0 we take two. for prelocate we take max.
eigFun = zeros(size(cosine,1),2*lastEigIdx); %if N==0 we take one fun, if N~=0 we take two. for prelocate we take max.
count = 1;
for i = 1:lastEigIdx
    order = C_idx(i)-1; % note we start from order zero;
    idxOfEig = R_idx(i);
    if gpu_use==1
        JrRhoGpu = gpuArray(JrRho(:,:,order+1)');
        JsampGpu = gpuArray(Jsamp(:,:,order+1));
        diagGpu = gpuArray(diag(rho.*psd.*quadKer));
        sqrtrSamprGpu = gpuArray(sqrtrSampr);
        eigFunTotGpu = eigFunTot(:,idxOfEig,order+1);
        quadNysGpu = quadNys;
        sqrtQuadNysGpu = gpuArray(sqrt(quadNys));
        Hsamp = sqrtrSamprGpu.*(JsampGpu*diagGpu*JrRhoGpu);
        vCorrect = (1./sqrtQuadNysGpu).*eigFunTotGpu;
        vNys = gather((Hsamp*(quadNysGpu.* vCorrect)) * (1/eigValTot(i)));
    else
        Hsamp = sqrtrSampr.*((Jsamp(:,:,order+1)*(diag(rho.*psd.*quadKer)*JrRho(:,:,order+1)')));
        vCorrect = (1./sqrt(quadNys)).*eigFunTot(:,idxOfEig,order+1); % correcting due to similarity. 
        vNys = (Hsamp*(quadNys.* vCorrect)) * (1/eigValTot(i)); 
    end
    if order==0
        vNys = reshape(vNys(idxRsamp),size(idxRsamp));
        eigFun(:,count) = (1/sqrt(2*pi))*vNys;
        eigVal(count) = eigValTot(i);
        count = count+1;
    else
        vNys = reshape(vNys(idxRsamp),size(idxRsamp));
        eigFun(:,count) = sqrt((1/pi))*vNys.*cosine(:,order+1);
        eigVal(count) = eigValTot(i);
        count = count+1;
        eigFun(:,count) = sqrt((1/pi))*vNys.*sine(:,order+1);
        eigVal(count) = eigValTot(i);
        count = count+1;
    end
end
 
eigVal = eigVal(eigVal>0); % we throw away the extra zeros;
eigFun = eigFun(:,1:length(eigVal)); % we throw away the extra zeros col;
end
    
    


    
        
        
       
