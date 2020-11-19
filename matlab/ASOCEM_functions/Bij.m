function b_ij = Bij(phi_ij,phi_ijp1,phi_im1j,phi_ip1j,mu,eta)
    b_ij = mu*(eta^2 + (phi_ijp1-phi_ij)^2 + ((phi_ip1j-phi_im1j)/2)^2)^(-0.5);
end