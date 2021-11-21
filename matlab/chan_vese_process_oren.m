function [phi,mu0_est,cov0_est,mu1_est,cov1_est] = chan_vese_process(I,phi_0,cov_mat_sz,dt,mu,nu,eta,Eps,maxIter,tol)
% This function compute the chan vese level set function phi under the assumption of radial autocorrelation.

  
%% time evolution process
% set boundary condition and start iteration
phi = NeumannBoundCondMod(phi_0);
for iter = 1:maxIter
    iter
    if nnz(phi>0)==0
        break
    end
    if nnz(phi<=0)==0
        break
    end
    % compute mean and cov
    % store patches
    patch_0 = [];
    patch_1 = [];
    area = sqrt(cov_mat_sz);
    for i=1:floor(size(I,1)/area)
        for j=1:floor(size(I,2)/area)
            tmp = phi((i-1)*area+1:i*area,(j-1)*area+1:j*area);           
            if nnz(tmp>=0)>0      
                patch_0 = [patch_0,reshape(I((i-1)*area+1:i*area,(j-1)*area+1:j*area),[],1)];
            end
            if nnz(tmp<=0)>0
                patch_1 = [patch_1,reshape(I((i-1)*area+1:i*area,(j-1)*area+1:j*area),[],1)];
            end
        end
    end
    if size(patch_0,2)<=10
        phi = ones(size(phi)); % we dont want to use this micrograph
        return
    end
    if size(patch_1,2)<=10
        phi = ones(size(phi)); % we dont want to use this micrograph
        return
    end
    % compte mean and cov
    mu0_est = mean(patch_0,2);
    cov0_est = cov(patch_0');
    mu1_est = mean(patch_1,2);
    cov1_est = cov(patch_1');
    % sizes derived from cov matrix
    cov0_inv = pinv(cov0_est);
    cov1_inv = pinv(cov1_est);
    logdet0 = logdetAmitay(cov0_est);
    logdet1 = logdetAmitay(cov1_est);
    area = sqrt(cov_mat_sz);
    % compute phi
    patch_r = 1;
    patch_c = 1;
    f = reshape(I(1:area,1:area),[],1);
    for i = 2 + floor(area/2) : size(phi,1)-1 - ceil(area/2)
        for j = 2 + floor(area/2) : size(phi,2)-1 - ceil(area/2) 
            if or(ceil(i/area) ~= patch_r,ceil(j/area) ~= patch_c)
                patch_r = ceil(i/area);
                patch_c = ceil(j/area);
                f = reshape(I((patch_r-1)*area+1:patch_r*area,(patch_c-1)*area+1:patch_c*area),[],1);      
            end
            phi(i,j) = timeEvCov(f,phi(i-1:i+1,j-1:j+1),mu,nu,...
              mu0_est,cov0_inv,logdet0,mu1_est,cov1_inv,logdet1,dt,eta,Eps);
        end
    end
    
    % boundary condition
    phi = NeumannBoundCondMod(phi); 
    
    % stopping criteria
    if mod(iter,5) == 0 && iter > 10
        area_new = phi>0; area_old = phi_old>0;
        changed_area =abs(area_new - area_old);
%         changedArea = (phi>0 | phiOld) - (phi>0 & phiOld); %union minus intersection
%         changedArea = sum(changedArea(:));

        if sum(changed_area(:))/sum(area_old(:)) < tol
%             fprintf('%%%%   stopped at %d-th iteration   %%%%\n',iter);
            break
        end
    end
    if mod(iter,5) == 0 
        phi_old = phi;
    end
end

end
