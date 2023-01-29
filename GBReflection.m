%% Reflection of the GB from the surface of a sphere (Frequency Independent)
function [Gamma_r_mat,eta_new_r1,eta_new_r2,sigma_r] =  GBReflection(xi_1,xi_2,nu_vect,eta_new_i1,eta_new_i2,sigma_i,Gamma_new_i_Q_Mat,l_b_vect)
% Step 4: Establishing a plane of incident and Rotation of the
% incident GB in the plane of incidence
theta_i_vect = acos(-dot(sigma_i,nu_vect,2)); % The initial angle of incidence; $$$$$$$$$ Chapter 3, equation (3.65) $$$$$$$$$
theta_i_mat = twodvect3dmat(theta_i_vect);
% xi_1, xi_2 and nu are the local coordinate on the sphere, at the point Q
% Establishing the matrix xi1xi2nu as a fucntion of the local
% coordinate system
Rot_lcs_xi_nu = [xi_1,xi_2,nu_vect];
Rot_lcs_xi_nu_mat_tr = twodvect3dmat(Rot_lcs_xi_nu);% $$$$$$$$$ Chapter 3, equation (3.60) $$$$$$$$$
%
% Matrix relation between (eta_i_1, eta_i_2, sigma_i) and (xi_1, xi_2, nu)
%Rot_etasig_nuxi = [dot(eta_new_i1,xi_1,2), dot(eta_new_i1,xi_2,2), dot(eta_new_i1,nu_vect,2),...
%dot(eta_new_i2,xi_1,2),dot(eta_new_i2,xi_2,2),dot(eta_new_i2,nu_vect,2),...
%dot(sigma_i,xi_1,2),dot(sigma_i,xi_2,2), dot(sigma_i,nu_vect,2)];
Rot_etasig_nuxi_11 = dot(eta_new_i1,xi_1,2);
Rot_etasig_nuxi_12 = dot(eta_new_i1,xi_2,2);
Rot_etasig_nuxi_13 = dot(eta_new_i1,nu_vect,2);
Rot_etasig_nuxi_1 = [Rot_etasig_nuxi_11,Rot_etasig_nuxi_12,Rot_etasig_nuxi_13];
%
Rot_etasig_nuxi_21 = dot(eta_new_i2,xi_1,2);
Rot_etasig_nuxi_22 = dot(eta_new_i2,xi_2,2);
Rot_etasig_nuxi_23 = dot(eta_new_i2,nu_vect,2);
Rot_etasig_nuxi_2 = [Rot_etasig_nuxi_21,Rot_etasig_nuxi_22,Rot_etasig_nuxi_23]; %Matrix relation between (eta_i_1,eta_i_2,sigma_i) and (xi_1,xi_2,nu)
%
Rot_etasig_nuxi_31 = dot(sigma_i,xi_1,2);
Rot_etasig_nuxi_32 = dot(sigma_i,xi_2,2); 
Rot_etasig_nuxi_33 = dot(sigma_i,nu_vect,2);
Rot_etasig_nuxi_3 = [Rot_etasig_nuxi_31,Rot_etasig_nuxi_32,Rot_etasig_nuxi_33];
Rot_etasig_nuxi = [Rot_etasig_nuxi_1,Rot_etasig_nuxi_2,Rot_etasig_nuxi_3];
Rot_etasig_i_nuxi_mat = twodvect3dmat(Rot_etasig_nuxi);
%
% Projection matrix of (eta_i_1,eta_i_2) on (xi_1,xi_2) $$$$$$$$$ Chapter 3, equation (3.25) and equation (3.26) $$$$$$$$$
Theta_eta_i = [Rot_etasig_nuxi_11, Rot_etasig_nuxi_12 ,Rot_etasig_nuxi_21, Rot_etasig_nuxi_22];
Theta_eta_i_vect = transpose(Theta_eta_i);
Theta_eta_i_vect_mat = pagetranspose(reshape(Theta_eta_i_vect,2,2,size(Theta_eta_i_vect,2)));
% Projection matrix of (eta_i_1,eta_i_2) on (xi_1,xi_2)
%Theta_eta_i = Rot_etasig_nuxi([1,2],[1,2]);
% Matrix relation between (eta_r_1, eta_r_2, sigma_i) and (xi_1, xi_2, nu) 
%Rot_etsig_nu_r = [-dot(eta_new_i1,xi_1,2), dot(eta_new_i1,xi_2,2), dot(eta_new_i1,nu_vect_inter,2);...
%dot(eta_new_i2,xi_1,2) dot(eta_new_i2,xi_2,2),dot(eta_new_i2,nu_vect_inter,2);...
%dot(sigma_i,xi_1,2),dot(sigma_i,xi_2,2), -dot(sigma_i,nu_vect_inter,2)]; 
Rot_etsig_nu_r_11 = -dot(eta_new_i1,xi_1,2);
Rot_etsig_nu_r_12 = dot(eta_new_i1,xi_2,2);
Rot_etsig_nu_r_13 = dot(eta_new_i1,nu_vect,2);
Rot_etasig_r_nuxi_1 = [Rot_etsig_nu_r_11,Rot_etsig_nu_r_12,Rot_etsig_nu_r_13];
%%%
Rot_etsig_nu_r_21 = dot(eta_new_i2,xi_1,2);
Rot_etsig_nu_r_22 = dot(eta_new_i2,xi_2,2);
Rot_etsig_nu_r_23 = dot(eta_new_i2,nu_vect,2);
Rot_etasig_r_nuxi_2 = [Rot_etsig_nu_r_21,Rot_etsig_nu_r_22,Rot_etsig_nu_r_23]; %Matrix relation between (eta_i_1,eta_i_2,sigma_i) and (xi_1,xi_2,nu)
%%%
Rot_etsig_nu_r_31 = dot(sigma_i,xi_1,2);
Rot_etsig_nu_r_32 = dot(sigma_i,xi_2,2);
Rot_etsig_nu_r_33 = -dot(sigma_i,nu_vect,2);
Rot_etasig_r_nuxi_3 = [Rot_etsig_nu_r_31,Rot_etsig_nu_r_32,Rot_etsig_nu_r_33];
Rot_etasig_r_nuxi = [Rot_etasig_r_nuxi_1,Rot_etasig_r_nuxi_2,Rot_etasig_r_nuxi_3];
Rot_etasig_r_nuxi_vect = transpose(Rot_etasig_r_nuxi);
Rot_etasig_r_nuxi_mat = pagetranspose(reshape(Rot_etasig_r_nuxi_vect,3,3,size(Rot_etasig_r_nuxi_vect,2)));
Rot_etasig_r_nuxi_mat = twodvect3dmat(Rot_etasig_r_nuxi);
% We used the relation sig_i.nu = - sig_r.nu, sig_i.xi2 = sig_r.xi2
% and etai2.nu = etar2.nu, etai2.xi2 = -etar2.xi2
%Reflected beam coordinates with respect to the global cordinates.
Rot_lbc_eta_sig_r_eta_mat = pagemtimes(Rot_etasig_r_nuxi_mat,Rot_lcs_xi_nu_mat_tr);%Matrix 
sigma_r = transpose(reshape(pagetranspose(Rot_lbc_eta_sig_r_eta_mat(3,:,:)),3,length(Rot_lbc_eta_sig_r_eta_mat)));
eta_new_r2 = transpose(reshape(pagetranspose(Rot_lbc_eta_sig_r_eta_mat(2,:,:)),3,length(Rot_lbc_eta_sig_r_eta_mat)));
eta_new_r1 = transpose(reshape(pagetranspose(Rot_lbc_eta_sig_r_eta_mat(1,:,:)),3,length(Rot_lbc_eta_sig_r_eta_mat)));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
J = [1 0;0 -1]; % Mirror symmetric matrixs
J_vect = [ones(length(l_b_vect),1),zeros(length(l_b_vect),1),zeros(length(l_b_vect),1),-ones(length(l_b_vect),1)];
J_mat = twodvect3dmat(J_vect);
c_radius = 1; 
c_11 = 1/c_radius; %These are only true for sphere. % One needs to change this for the different geometry Ofcourse. 
c_12 = 0;
c_21 = 0;
c_22 = 1/c_radius;
C = [c_11,c_12,c_21,c_22]; % this is only true for a sphere, for a general object it needs to be calculate locally
C_vect = [c_11 + zeros(length(l_b_vect),1),c_12 + zeros(length(l_b_vect),1),c_21 + zeros(length(l_b_vect),1),c_22 + zeros(length(l_b_vect),1)];
C_mat = twodvect3dmat(C_vect);
% Calculating the inverse of the Theta_Matrix
Theta_eta_i = [Rot_etasig_nuxi_11, Rot_etasig_nuxi_12 ,Rot_etasig_nuxi_21, Rot_etasig_nuxi_22];
det_Theta_eta_i = Theta_eta_i(:,1).*Theta_eta_i(:,4) - Theta_eta_i(:,2).*Theta_eta_i(:,3);
inv_Theta_eta_i = [Theta_eta_i(:,4),-Theta_eta_i(:,2),-Theta_eta_i(:,3),Theta_eta_i(:,1)];
inv_Theta_eta_i = inv_Theta_eta_i./det_Theta_eta_i; % Inverse of Theta Matrix
inv_Theta_eta_i_mat = pageinv(Theta_eta_i_vect_mat);
%
%Calculating the Reflecetd Gaussian beam strucutre
theta_J_mat = pagemtimes(inv_Theta_eta_i_mat,J_mat); 
% theta_J_mat_trans = pagetranspose(theta_J_mat);
Gamma_r_mat = pagemtimes(J_mat, pagemtimes(Gamma_new_i_Q_Mat,J_mat))  + 2*cos(theta_i_mat).*pagemtimes(theta_J_mat,pagemtimes(C_mat,theta_J_mat));
% 
sig_r_Q = sigma_i - 2*(dot(sigma_i,nu_vect,2)).*nu_vect; %Computing the reflected field direction, the reflected field shall acts as the new incident GB for our algorithm
sig_r_Q = sig_r_Q./vecnorm(sig_r_Q,2,2);