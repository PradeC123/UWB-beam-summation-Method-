function [Gamma_new_i_Q_Mat,A_i_Q,xi_1,xi_2,eta_new_i1,eta_new_i2,sigma_i] =  Gamma_i_rot(Z,l_b_vect,b,~,~,nu_vect)
%% Estabilishing the Plane of incidence and the rotation of the incident GB in the plane of incidence(Frequency Independent)
% Gamma_i and its inverse. 
theta_n = Z(:,6);
phi_n = Z(:,7);
sigma_i = [sin(theta_n).*cos(phi_n),sin(theta_n).*sin(phi_n),cos(theta_n)]; %New beams in the sigma direction 
eta_i1 = [cos(theta_n).*cos(phi_n),cos(theta_n).*sin(phi_n),-sin(theta_n)]; %transvere direction of the GB 
eta_i2 = [-sin(phi_n),cos(phi_n),zeros(length(phi_n),1)]; %transverse direction of GB 
Gamma_i = [1./(1j*b*cos(theta_n).^2),zeros(length(theta_n),1),zeros(length(theta_n),1),zeros(length(theta_n),1)+1/(1j*b)]; 
%Complex curvature matrix of GB. 
det_Gamma_i = Gamma_i(:,1).*Gamma_i(:,4) - Gamma_i(:,2).*Gamma_i(:,3);
inv_Gamma_i = [Gamma_i(:,4),-Gamma_i(:,2),-Gamma_i(:,3),Gamma_i(:,1)];
inv_Gamma_i = inv_Gamma_i./det_Gamma_i; 
Gamma_i_tr_1 = twodvect3dmat(Gamma_i);
inv_Gamma_i_mat_1 = pageinv(Gamma_i_tr_1); %Inverse of Gamma_i
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Gamma_Q, A_i_Q at the specular point $$$$$$$$$ Chapter 3, equation (3.2) $$$$$$$$$
l_b_matrix = [l_b_vect,zeros(length(l_b_vect),1),zeros(length(l_b_vect),1),l_b_vect];
inv_Gamma_Q =  inv_Gamma_i + l_b_matrix; 
L_b_vect = [l_b_vect,zeros(length(l_b_vect),1),zeros(length(l_b_vect),1),l_b_vect];
L_b_mat_tr = twodvect3dmat(L_b_vect);
inv_Gamma_i_Mat = inv_Gamma_i_mat_1+L_b_mat_tr;
Gamma_i_Q_Mat = pageinv(inv_Gamma_i_Mat);
det_inv_Gamma_Q = inv_Gamma_Q(:,1).*inv_Gamma_Q(:,4) - inv_Gamma_Q(:,2).*inv_Gamma_Q(:,3);
Gamma_i_Q = [inv_Gamma_Q(:,4),-inv_Gamma_Q(:,2),-inv_Gamma_Q(:,3),inv_Gamma_Q(:,1)];
Gamma_i_Q = Gamma_i_Q./det_inv_Gamma_Q; % Complex Curvatuer matrix at the specular point Q
det_Gamma_Q = Gamma_i_Q(:,1).*Gamma_i_Q(:,4) - Gamma_i_Q(:,2).*Gamma_i_Q(:,3); % determinant ot the Gamma_Q matrix at the specular point 
A_i_Q = sqrt(det_Gamma_Q./det_Gamma_i); % Complex Amplitude of the Beam at the specular point
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Establishing the Local Plane of incidence $$$$$$$$$ Chapter 3, equation (3.67) $$$$$$$$$
xi_2 = cross(nu_vect,-sigma_i,2); % The vector perpendicular to the plane of incidence 
indx_xi_2_0 = find(vecnorm(xi_2,2,2)==0);
xi_2(indx_xi_2_0,:) = eta_i2(indx_xi_2_0,:); % This case implies Plane of incidence doesnot exists and nu_vect and sigma are coliner
xi_2 = xi_2./vecnorm(xi_2,2,2); %Normalization of the vector perpendicular to the plane of incidence
xi_1 = cross(xi_2,nu_vect,2);
xi_1 = xi_1./vecnorm(xi_1,2,2); %Normalization of the vector in the plane of incidence
%Choosing the desired beam coordinate system such that, eta_new_i1 is in
% the perpendicular direction of the plane of incidence and eta_new_i2 is in the plane of
% incendence
eta_new_i2 = xi_2;
eta_new_i1 = cross(eta_new_i2,sigma_i,2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The matrix that rotates the plane of the transverse coordinate system to
% the new desired plane such that it satisfies the specified conditions as
% mentioned above. $$$$$$$$$ Chapter 3, equation (3.69) $$$$$$$$$
theta_eta_12 = [ dot(eta_new_i1,eta_i1,2), dot(eta_new_i1,eta_i2,2), dot(eta_new_i2,eta_i1,2), dot(eta_new_i2,eta_i2,2)];
theta_eta_12_mat_tr = twodvect3dmat(theta_eta_12);
% $$$$$$$$$ Chapter 3, equation (3.71) $$$$$$$$$
Gamma_new_i_Q_Mat = pagemtimes(theta_eta_12_mat_tr,pagemtimes(Gamma_i_Q_Mat,pagetranspose(theta_eta_12_mat_tr)));%Matrix
%multiplication 
end