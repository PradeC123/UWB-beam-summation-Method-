%% Rotation of the transverse coordinate (Frequency Independent)
function [Gamma_new_r_mat,Gamma_new_r,theta_sig_r,phi_sig_r] = N2FF_Gammar_rot(Gamma_r_mat,eta_new_r1,eta_new_r2,sigma_r)
%Step 5: Rotation of the transverse coordinate 
theta_sig_r = acos(sigma_r(:,3));
phi_sig_r = mod(atan2(sigma_r(:,2),sigma_r(:,1)),2*pi);
eta_tilde_1r = [cos(theta_sig_r).*cos(phi_sig_r),cos(theta_sig_r).*sin(phi_sig_r),-sin(theta_sig_r)];
eta_tilde_2r = [-sin(phi_sig_r),cos(phi_sig_r), zeros(length(phi_sig_r),1)];
Phi_eta = [dot(eta_tilde_1r,eta_new_r1,2),dot(eta_tilde_1r,eta_new_r2,2),dot(eta_tilde_2r,eta_new_r1,2),dot(eta_tilde_2r,eta_new_r2,2)];
Phi_eta_mat = twodvect3dmat(Phi_eta);
% Roation of the Gamma_i_Q in the plane of incidence
Gamma_new_r_mat = pagemtimes(Phi_eta_mat,pagemtimes(Gamma_r_mat,pagetranspose(Phi_eta_mat)));
Gamma_new_r = [transpose(reshape(pagetranspose(Gamma_new_r_mat(1,:,:)),2,length(Gamma_new_r_mat))), ...
    transpose(reshape(pagetranspose(Gamma_new_r_mat(2,:,:)),2,length(Gamma_new_r_mat)))]; 
end