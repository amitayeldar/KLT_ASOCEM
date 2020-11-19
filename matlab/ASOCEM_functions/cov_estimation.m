function [cov_mat,R_est] = cov_estimation(I,phi,idx_mat,max_d,mean_region,sign)


% set region of interest
if sign == 1
    region = phi > 0;
end
if sign == -1
    region = phi <= 0; 
end
if nnz(region)<1
    error('phi converged to trivial case')
end

% make I=image "square"
if size(I,1)~=size(I,2)
   tmp1 = size(I,1);
   tmp2 = size(I,2);
   if size(I,1)>size(I,2)
       I(tmp2+1:tmp1,tmp2+1:tmp1) = 0;
       region(tmp2+1:tmp1,tmp2+1:tmp1) = 0;
   else
       I(tmp1+1:tmp2,tmp1+1:tmp2) = 0;
       region(tmp1+1:tmp2,tmp1+1:tmp2) = 0;
   end
end


 samples_idx = find(region>0); % estimate autocorr for pixles in region only
 I(region>0) = I(region>0) - mean_region; % recentre mean to zero;
 
 [R_est,~,cnt]=cryo_epsdR(I,samples_idx,max_d);
 if size(R_est,1)<max(idx_mat(:))
     error('region is not big enough to estimate cov');
 end
    cov_mat = R_est(idx_mat(:));
    cov_mat = reshape(cov_mat,sqrt(size(cov_mat,1)),sqrt(size(cov_mat,1)));

    

end

