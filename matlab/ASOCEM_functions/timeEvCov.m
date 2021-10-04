function phi_ij = timeEvCov(f,PHI,mu,nu,mu0_est,cov0_inv,logdet0,mu1_est,cov1_inv,logdet1,dt,eta,Eps)
%% This function is used for the time evaluation of phi as written in the submited paper
% compute terms for the time evaluation
c = 2;
a_im1j = Aij(PHI(c-1,c),PHI(c,c),PHI(c-1,c-1),PHI(c-1,c+1),mu,eta);
a_ij = Aij(PHI(c,c),PHI(c+1,c),PHI(c,c-1),PHI(c,c+1),mu,eta);
b_ijm1 = Bij(PHI(c,c-1),PHI(c,c),PHI(c-1,c-1),PHI(c+1,c-1),mu,eta);
b_ij = Bij(PHI(c,c),PHI(c,c+1),PHI(c-1,c),PHI(c+1,c),mu,eta);
delta_ij = deltaEps(PHI(c,c),Eps);

% time evaluation
% fidelity terms
RT = - nu +(1/size(f,1))*0.5*((logdet1-logdet0)+(f-mu1_est)'*cov1_inv*(f-mu1_est)-(f-mu0_est)'*cov0_inv*(f-mu0_est));
% curvature flow terms
CFT = a_ij*PHI(c+1,c) + a_im1j*PHI(c-1,c) + b_ij*PHI(c,c+1) + b_ijm1*PHI(c,c-1); 

phi_ij = ( PHI(c,c) + dt * delta_ij * (CFT + RT ) )...
		 /( 1 + dt * delta_ij * ( a_ij + a_im1j + b_ij + b_ijm1 ) );
		 
     
end
