function [phi,mu0_est,cov0_est,mu1_est,cov1_est] = chan_vese_process(I,phi_0,cov_mat_sz,dt,mu,nu,eta,Eps,maxIter,tol,fast_flag)
% This function compute the chan vese level set function phi under the assumption of radial autocorrelation.

%% time evolution process
% set boundary condition and start iteration
phi = NeumannBoundCondMod(phi_0);
for iter = 1:maxIter
    if nnz(phi>0)==0
        break
    end
    if nnz(phi<=0)==0
        break
    end
    % compute mean and cov
    
    % store patches
    cnt_0 = 1;
    cnt_1 = 1;

    area = sqrt(cov_mat_sz);
    patch_0 = zeros(area^2,1);
    patch_1 = zeros(area^2,1);
    for i=1:floor(size(I,1)/area)
        for j=1:floor(size(I,2)/area)
            tmp = phi((i-1)*area+1:i*area,(j-1)*area+1:j*area);           
            if nnz(tmp>=0)>0   
                patch_0(:,cnt_0) = reshape(I((i-1)*area+1:i*area,(j-1)*area+1:j*area),[],1);
                cnt_0 = cnt_0 +1;
            end
            if nnz(tmp<=0)>0
                patch_1(:,cnt_1) = reshape(I((i-1)*area+1:i*area,(j-1)*area+1:j*area),[],1);
                cnt_1 = cnt_1 +1;
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
    %% compute phi
    if fast_flag==0
        % compute phi near zero
        band_width = min(min(min(phi(phi>0))),abs(max(max(phi(phi<0)))));
        stop = 0;
        while stop ==0
            tmp_p = and(phi>=0,phi<band_width);
            tmp_m = and(phi<0,phi>-band_width);
            if nnz(tmp_p)+nnz(tmp_m)>0.5*size(phi,1)*size(phi,2)
                break
            end
            band_width = 2*band_width;
        end
        [row,col] = find(and(phi<band_width,phi>-band_width)); 
        for i=1:size(row,1)
            if and(and(2<=row(i),row(i)<=size(phi,1)-rem(size(phi,1),area)-1),and(2<=col(i),col(i)<=size(phi,2)-rem(size(phi,2),area)-1))==1
                    patch_r = ceil(row(i)/area);
                    patch_c = ceil(col(i)/area);
                    f = reshape(I((patch_r-1)*area+1:patch_r*area,(patch_c-1)*area+1:patch_c*area),[],1);      
                    phi(row(i),col(i)) = timeEvCov(f,phi(row(i)-1:row(i)+1,col(i)-1:col(i)+1),mu,nu,...
                               mu0_est,cov0_inv,logdet0,mu1_est,cov1_inv,logdet1,dt,eta,Eps);
            end
        end
    end
    if fast_flag==1
        phi_old = phi;
        for i=2:size(phi,1)-rem(size(phi,1),area)-1
            parfor j=2:size(phi,2)-rem(size(phi,2),area)-1
                patch_r = ceil(i/area);
                patch_c = ceil(j/area);
                f = reshape(I((patch_r-1)*area+1:patch_r*area,(patch_c-1)*area+1:patch_c*area),[],1);                  
                phi(i,j)=timeEvCov(f,phi_old(i-1:i+1,j-1:j+1),mu,nu,...
                               mu0_est,cov0_inv,logdet0,mu1_est,cov1_inv,logdet1,dt,eta,Eps);
            end
        end
    end
    
    if fast_flag==2
        for i = 1:area:size(phi,1)-rem(size(phi,1),area)-1
            for j  = 1:area:size(phi,2)-rem(size(phi,2),area)-1
                    f = reshape(I(i:area+i-1,j:area+j-1),[],1); 
                    delta_ij = deltaEps(phi(i,j),Eps);
                    RT = - nu +(1/size(f,1))*0.5*((logdet1-logdet0)+(f-mu1_est)'*cov1_inv*(f-mu1_est)-(f-mu0_est)'*cov0_inv*(f-mu0_est));
                    phi(i:area+i-1,j:area+j-1) = phi(i:area+i-1,j:area+j-1) + dt * delta_ij * RT; 
%                     phi(i:area+i-1,j:area+j-1) =  sign(RT);
            end
        end
    end
    % boundary condition
    phi = NeumannBoundCondMod(phi); 
    
    % stopping criteria
    if mod(iter,5) == 0 && iter > 10
        area_new = phi>0; area_old = phi_old_stop>0;
        changed_area =abs(area_new - area_old);
        if sum(changed_area(:))/sum(area_old(:)) < tol
            break
        end
    end
    if mod(iter,5) == 0 
        phi_old_stop = phi;
    end
end

end
