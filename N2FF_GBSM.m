%% Monostatic RCS Calculations in the far zone.(Frequency Dependent)
function [Sph_mnst,Z_r,N] = N2FF_GBSM(Gamma_r_mat,Gamma_new_r,l_b_vect,theta_obser,phi_obser,phi_sig_r,theta_sig_r,Q_S,k,A_i_Q,Reconcoeff,Z_r)
r_hat = [zeros(length(phi_sig_r),1) + sin(theta_obser)*cos(phi_obser),...
    zeros(length(phi_sig_r),1) + sin(theta_obser)*sin(phi_obser),...
    zeros(length(phi_sig_r),1) + cos(theta_obser)]; % Monostatic observation angle
theta = acos(r_hat(:,3));
r = 10^8;
phi = mod(atan2(r_hat(:,2),r_hat(:,1)),2*pi);
sig_dir = sin(theta_sig_r).*cos(phi_sig_r).*(r*sin(theta).*cos(phi)-Q_S(:,1)) + sin(theta_sig_r).*sin(phi_sig_r).*(r*sin(theta).*sin(phi)-Q_S(:,2))...
         + cos(theta_sig_r).*(r*cos(theta)-Q_S(:,3));
a_11 = sin(theta_obser).*cos(theta_sig_r).*cos(phi_obser-phi_sig_r)-cos(theta_obser).*sin(theta_sig_r);
b_11 = sin(theta_obser).*sin(phi_obser-phi_sig_r);
thet_vect = [a_11,zeros(length(a_11),1),zeros(length(l_b_vect),1),b_11];
thet_nat = twodvect3dmat(thet_vect);%
%I_vect = [ones(length(a_11),1),zeros(length(a_11),1),zeros(length(l_b_vect),1),ones(length(a_11),1)];
%I_nat = twodvect3dmat(I_vect);
%
det_Gamma_new_r = Gamma_new_r(:,1).*Gamma_new_r(:,4) - Gamma_new_r(:,2).*Gamma_new_r(:,3);
inv_Gamma_new_r = [Gamma_new_r(:,4),-Gamma_new_r(:,2),-Gamma_new_r(:,3),Gamma_new_r(:,1)];
inv_Gamma_new_r = inv_Gamma_new_r./det_Gamma_new_r; 
inv_Gamma_new_r_mat = pageinv(Gamma_r_mat);
exp_ff_GB = sum(sum(pagemtimes(thet_nat,pagemtimes(inv_Gamma_new_r_mat,pagetranspose(thet_nat)))));
exp_ff_GB = reshape(exp_ff_GB,length(exp_ff_GB),1);
freq_ind_par = exp_ff_GB/2 + dot(r_hat,Q_S,2) - l_b_vect;
if ~isempty(Z_r) 
  Z_r(:,length(Z_r(1,:))+1) = freq_ind_par;
  Z_r(:,length(Z_r(1,:))+1) = A_i_Q.*sqrt(1./(det_Gamma_new_r)).*((sig_dir>0.85*r)==1);
  %Z_r(:,length(Z_r(1,:))+1) = Reconcoeff.*A_i_Q.*sqrt(1./(det_Gamma_new_r)).*exp(1j*k*freq_ind_par).*((sig_dir>0.85*r)==1);
end
Sph_mnst =  sum(Reconcoeff.*A_i_Q.*sqrt(1./(det_Gamma_new_r)).*exp(1j*k*freq_ind_par).*((sig_dir>0.9*r)==1));
end