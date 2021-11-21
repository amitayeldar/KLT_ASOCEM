function a_ij = Aij(phi_ij,phi_ip1j,phi_ijm1,phi_ijp1,mu,eta)
    a_ij = mu*(eta^2 + (phi_ip1j-phi_ij)^2 + ((phi_ijp1-phi_ijm1)/2)^2)^(-0.5);
end