function [phi,mu0,R0_est,mu1,R1_est] = chan_vese_process(I,phi_0,cov_mat_sz,dt,mu,nu,eta,Eps,maxIter,tol)

%% initial steps

max_d = floor(sqrt(2)*sqrt(cov_mat_sz)) + 1;

%% from 1d isotropic to covarience matrix
% first arrange an R^2 index matrix
idx_cell= cell(cov_mat_sz,1);
cnt = 1;
for j=1:sqrt(cov_mat_sz)
    for i=1:sqrt(cov_mat_sz)
       idx_cell{cnt} = [i,j];
       cnt = cnt + 1;
    end
end
% second construct index matrix for cov_rad_mat
mat_for_estimate_r_cov_samp_point = ones(max_d,max_d);
[~,x,~]=cryo_epsdR(mat_for_estimate_r_cov_samp_point,1:max_d^2,max_d);
idx_mat = zeros(cov_mat_sz,cov_mat_sz);
for i = 1:cov_mat_sz
    for j = 1:cov_mat_sz
        d_ij = norm(idx_cell{i}-idx_cell{j});
        idx_mat(i,j) = find(x==d_ij);
    end
end

% set boundary condition and start iteration
phi = NeumannBoundCondMod(phi_0);  



%% evolution process
for iter = 1:maxIter
    iter
    % compute expectation
    if nnz(phi>0)==0
        break
    end
    if nnz(phi<=0)==0
        break
    end
    mu0 = mean( I( phi > 0 ) );
    mu1 = mean( I( phi <= 0 ) );
    
    %compute covariance
    [s0,R0_est] = cov_estimation(I,phi,idx_mat,max_d,mu0,1);
    [s1,R1_est] = cov_estimation(I,phi,idx_mat,max_d,mu1,-1);
    % sizes derived from cov matrix
    s1_inv = pinv(s0);
    s2_inv = pinv(s1);
    logdet0 = logdetAmitay(s0);
    logdet1 = logdetAmitay(s1);
    area = sqrt(cov_mat_sz);
    % compute phi
    for i = 2 + floor(area/2) : size(phi,1)-1 - floor(area/2)
        for j = 2 + floor(area/2) : size(phi,2)-1 - floor(area/2)  
            f = I(i-floor(area/2):i+floor(area/2),j-floor(area/2):j+floor(area/2));        
            f = f(:);
            phi(i,j) = timeEvCov(f,phi(i-1:i+1,j-1:j+1),mu,nu,...
              mu0,s1_inv,logdet0,mu1,s2_inv,logdet1,dt,eta,Eps);
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
            fprintf('%%%%   stopped at %d-th iteration   %%%%\n',iter);
            break
        end
    end
    if mod(iter,5) == 0 
        phi_old = phi;
    end
end

end
