%% Expansion paramters
a = 1;
b = 1;
k_max = 1000;
nu_max = 0.20;
z_explane = 0;
k = 1000;
theta_inc = 0;
phi_inc = 0;
%% Parameters of the sphere
c_1 = [0,0,0]; %The location of the center of the sphere 
r_1 = 1; % The radius of the sphere 
%% Establishing the Expansion Plane (GB)
tic
Expplane = Chapter4Expansionplane(a,b,k_max,nu_max,z_explane);
Z = Expplane;
toc
%% Calculations of Reconstruction Coeff
tic
[Reconcoeff,Z] = Chapter4Reconstructioncoeff(Z,theta_inc,phi_inc,nu_max,k_max,k,b,z_explane);
toc
%% Intersection of a beam lattice with a Sphere 
tic
theta_n = Z(:,6);
phi_n = Z(:,7);
x_bar = sqrt(2*pi*nu_max*b/k_max);
sigma_i = [sin(theta_n).*cos(phi_n),sin(theta_n).*sin(phi_n),cos(theta_n)];
X_0 = [Z(:,3)*x_bar,Z(:,4)*x_bar, zeros(length(Z(:,3)),1)+z_explane];
Delta_1 = dot(sigma_i,X_0-c_1,2).^2 - vecnorm(X_0-c_1,2,2).^2 + r_1^2;
Z(find(Delta_1<0),:) = []; %Removing the beams that doesnot hit the sphere
Reconcoeff(find(Delta_1<0),:) = [];
Delta_1(find(Delta_1<0),:) = []; %Removing the beams that doesnot hit the sphere
% Post filteration process 
theta_n = Z(:,6);
phi_n = Z(:,7);
sigma_i = [sin(theta_n).*cos(phi_n),sin(theta_n).*sin(phi_n),cos(theta_n)]; %New beams in the sigma direction 
eta_i1 = [cos(theta_n).*cos(phi_n),cos(theta_n).*sin(phi_n),-sin(theta_n)]; %transvere direction of the GB 
eta_i2 = [-sin(phi_n),cos(phi_n),zeros(length(phi_n),1)]; %transverse direction of GB 
X_0 = [Z(:,3)*x_bar,Z(:,4)*x_bar, zeros(length(Z(:,3)),1)+z_explane]; %New beams in the X_0 expansion plane direction 
d_1 = -dot(sigma_i,X_0-c_1,2) + sqrt(Delta_1);
d_2 = -dot(sigma_i,X_0-c_1,2) - sqrt(Delta_1);
d = min(d_1,d_2); 
Q_S = X_0 + d.*sigma_i; % Specular point on the sphere
nu_vect = (Q_S-c_1)./vecnorm(Q_S-c_1,2,2); %Local Outward vector
l_b_vect = vecnorm(Q_S-X_0,2,2); %Distance from the expansion plane to the specular point 
%% Estabilishing the Plane of incidence and the rotation of the incident GB in the plane of incidence
% Gamma_i and its inverse. 
Gamma_i = [1./(1j*b*cos(theta_n).^2),zeros(length(theta_n),1),zeros(length(theta_n),1),zeros(length(theta_n),1)+1/(1j*b)]; %Complex curvature matrix of GB. 
det_Gamma_i = Gamma_i(:,1).*Gamma_i(:,4) - Gamma_i(:,2).*Gamma_i(:,3);
inv_Gamma_i = [Gamma_i(:,4),-Gamma_i(:,2),-Gamma_i(:,3),Gamma_i(:,1)];
inv_Gamma_i = inv_Gamma_i./det_Gamma_i; 
%Gamma_i_tr = transpose(Gamma_i);
%Gamma_i_tr_1 = pagetranspose(reshape(Gamma_i_tr,2,2,size(Gamma_i,2)));
%inv_Gamma_i_mat_1 = pageinv(Gamma_i_tr_1); %Inverse of Gamma_i
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Gamma_Q, A_i_Q at the specular point
l_b_matrix = [l_b_vect,zeros(length(l_b_vect),1),zeros(length(l_b_vect),1),l_b_vect];
inv_Gamma_Q =  inv_Gamma_i + l_b_matrix; 
%L_b_vect = [l_b_vect,zeros(length(l_b_vect),1),zeros(length(l_b_vect),1),l_b_vect];
%L_b_vect_tr = transpose(L_b_vect);
%L_b_mat_tr = pagetranspose(reshape(L_b_vect_tr,2,2,size(L_b_vect_tr,2)));
%inv_Gamma_i_Mat = inv_Gamma_i_mat_1+L_b_mat_tr;
%Gamma_i_Q_Mat = pageinv(inv_Gamma_i_Mat);
det_inv_Gamma_Q = inv_Gamma_Q(:,1).*inv_Gamma_Q(:,4) - inv_Gamma_Q(:,2).*inv_Gamma_Q(:,3);
Gamma_i_Q = [inv_Gamma_Q(:,4),-inv_Gamma_Q(:,2),-inv_Gamma_Q(:,3),inv_Gamma_Q(:,1)];
Gamma_i_Q = Gamma_i_Q./det_inv_Gamma_Q; % Complex Curvatuer matrix at the specular point Q
det_Gamma_Q = Gamma_i_Q(:,1).*Gamma_i_Q(:,4) - Gamma_i_Q(:,2).*Gamma_i_Q(:,3); % determinant ot the Gamma_Q matrix at the specular point 
A_i_Q = sqrt(det_Gamma_Q./det_Gamma_i); % Complex Amplitude of the Beam at the specular point
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Establishing the Local Plane of incidence
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
% mentioned above.
theta_eta_12 = [ dot(eta_new_i1,eta_i1,2), dot(eta_new_i1,eta_i2,2), dot(eta_new_i2,eta_i1,2), dot(eta_new_i2,eta_i2,2)];
a_11 = dot(eta_new_i1,eta_i1,2); 
b_11 = dot(eta_new_i1,eta_i2,2);
c_11 = dot(eta_new_i2,eta_i1,2);
d_11 = dot(eta_new_i2,eta_i2,2);
Gamma_i_Q_1 = Gamma_i_Q(:,1); %Gamma_i_11
Gamma_i_Q_2 = Gamma_i_Q(:,2); %Gamma_i_12
Gamma_i_Q_3 = Gamma_i_Q(:,3); %Gamma_i_21
Gamma_i_Q_4 = Gamma_i_Q(:,4); %Gamma_i_22
%theta_eta_12_vect_tr = transpose(theta_eta_12);
%theta_eta_12_mat_tr = pagetranspose(reshape(theta_eta_12_vect_tr,2,2,size(theta_eta_12_vect_tr,2)));
%Gamma_new_i_Q_Mat = pagemtimes(theta_eta_12_mat_tr,pagemtimes(Gamma_i_Q_Mat,pagetranspose(theta_eta_12_mat_tr)));%Matrix
%multiplication 
% Roation of the Gamma_i_Q in the plane of incidence
Gamma_new_i_Q_1 = a_11.^2.*Gamma_i_Q_1 + a_11.*b_11.*Gamma_i_Q_2 + a_11.*b_11.*Gamma_i_Q_3 + b_11.^2.*Gamma_i_Q_4;
Gamma_new_i_Q_2 = a_11.*c_11.*Gamma_i_Q_1 + a_11.*d_11.*Gamma_i_Q_2 + c_11.*b_11.*Gamma_i_Q_3 + b_11.*d_11.*Gamma_i_Q_4;
Gamma_new_i_Q_3 = a_11.*c_11.*Gamma_i_Q_1 + c_11.*b_11.*Gamma_i_Q_2 + a_11.*d_11.*Gamma_i_Q_3 + b_11.*d_11.*Gamma_i_Q_4;
Gamma_new_i_Q_4 = c_11.^2.*Gamma_i_Q_1 + c_11.*d_11.*Gamma_i_Q_2 + c_11.*d_11.*Gamma_i_Q_3 + d_11.^2.*Gamma_i_Q_4;
Gamma_new_i_Q = [Gamma_new_i_Q_1,Gamma_new_i_Q_2,Gamma_new_i_Q_3,Gamma_new_i_Q_4];% Rotation of the Gamma_i_Q in the plane of incidence
%% Reflection of the GB from the surface of a sphere
% Step 4: Establishing a plane of incident and Rotation of the
% incident GB in the plane of incidence
theta_i_vect = acos(-dot(sigma_i,nu_vect,2)); % The initial angle of incidence
% xi_1, xi_2 and nu are the local coordinate on the sphere, at the point Q
% Establishing the matrix xi1xi2nu as a fucntion of the local
% coordinate system
Rot_lcs_xi_nu_11 = xi_1(:,1);
Rot_lcs_xi_nu_12 = xi_1(:,2);
Rot_lcs_xi_nu_13 = xi_1(:,3);
%
Rot_lcs_xi_nu_21 = xi_2(:,1);
Rot_lcs_xi_nu_22 = xi_2(:,2);
Rot_lcs_xi_nu_23 = xi_2(:,3);
%
Rot_lcs_xi_nu_31 = nu_vect(:,1); 
Rot_lcs_xi_nu_32 = nu_vect(:,2); 
Rot_lcs_xi_nu_33 = nu_vect(:,3); 
%Rot_lcs_xi_nu = [xi_1, xi_2,nu_vect];
%Rot_lcs_xi_nu_vect = transpose(Rot_lcs_xi_nu);
%Rot_lcs_xi_nu_mat = pagetranspose(reshape(Rot_lcs_xi_nu_vect,3,3,size(Rot_lcs_xi_nu_vect,2)));
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
%Rot_etasig_nuxi_vect = transpose(Rot_etasig_nuxi);
%Rot_etasig_nuxi_mat = pagetranspose(reshape(Rot_etasig_nuxi_vect,3,3,size(Rot_etasig_nuxi_vect,2)));
%
% Projection matrix of (eta_i_1,eta_i_2) on (xi_1,xi_2) 
Theta_eta_i = [Rot_etasig_nuxi_11, Rot_etasig_nuxi_12 ,Rot_etasig_nuxi_21, Rot_etasig_nuxi_22];
%Theta_eta_i_vect = transpose(Theta_eta_i);
%Theta_eta_i_vect_mat = pagetranspose(reshape(Theta_eta_i_vect,2,2,size(Theta_eta_i_vect,2)));
% Projection matrix of (eta_i_1,eta_i_2) on (xi_1,xi_2)
%Theta_eta_i = Rot_etasig_nuxi([1,2],[1,2]);
% Matrix relation between (eta_r_1, eta_r_2, sigma_i) and (xi_1, xi_2, nu) 
%Rot_etsig_nu_r = [-dot(eta_new_i1,xi_1,2), dot(eta_new_i1,xi_2,2), dot(eta_new_i1,nu_vect_inter,2);...
%dot(eta_new_i2,xi_1,2) dot(eta_new_i2,xi_2,2),dot(eta_new_i2,nu_vect_inter,2);...
%dot(sigma_i,xi_1,2),dot(sigma_i,xi_2,2), -dot(sigma_i,nu_vect_inter,2)]; 
Rot_etsig_nu_r_11 = -dot(eta_new_i1,xi_1,2);
Rot_etsig_nu_r_12 = dot(eta_new_i1,xi_2,2);
Rot_etsig_nu_r_13 = dot(eta_new_i1,nu_vect,2);
%
Rot_etsig_nu_r_21 = dot(eta_new_i2,xi_1,2);
Rot_etsig_nu_r_22 = dot(eta_new_i2,xi_2,2);
Rot_etsig_nu_r_23 = dot(eta_new_i2,nu_vect,2);
% 
Rot_etsig_nu_r_31 = dot(sigma_i,xi_1,2);
Rot_etsig_nu_r_32 = dot(sigma_i,xi_2,2);
Rot_etsig_nu_r_33 = -dot(sigma_i,nu_vect,2);
% We used the relation sig_i.nu = - sig_r.nu, sig_i.xi2 = sig_r.xi2
% and etai2.nu = etar2.nu, etai2.xi2 = -etar2.xi2
%Reflected beam coordinates with respect to the global cordinates.
%Rot_lbc_eta_sig_r = Rot_etsig_nu_r * Rot_lcs_xi_nu;
%Rot_lbc_eta_sig_r_mat = pagetranspose(reshape(Rot_etasig_nuxi_vect,3,3,size(Rot_etasig_nuxi_vect,2)));
%pagemtimes(Rot_etasig_nuxi_mat,Rot_lcs_xi_nu_mat);%Matrix %Continue from
%here tomorrow
%
Rot_lbc_eta_sig_r_eta_r1_11 = Rot_etsig_nu_r_11.*Rot_lcs_xi_nu_11 + Rot_etsig_nu_r_12.*Rot_lcs_xi_nu_21 + Rot_etsig_nu_r_13.*Rot_lcs_xi_nu_31;
Rot_lbc_eta_sig_r_eta_r1_12 = Rot_etsig_nu_r_11.*Rot_lcs_xi_nu_12 + Rot_etsig_nu_r_12.*Rot_lcs_xi_nu_22 + Rot_etsig_nu_r_13.*Rot_lcs_xi_nu_32;
Rot_lbc_eta_sig_r_eta_r1_13 = Rot_etsig_nu_r_11.*Rot_lcs_xi_nu_13 + Rot_etsig_nu_r_12.*Rot_lcs_xi_nu_23 + Rot_etsig_nu_r_13.*Rot_lcs_xi_nu_33;
%
Rot_lbc_eta_sig_r_eta_r2_21 = Rot_etsig_nu_r_21.*Rot_lcs_xi_nu_11 + Rot_etsig_nu_r_22.*Rot_lcs_xi_nu_21 + Rot_etsig_nu_r_23.*Rot_lcs_xi_nu_31;
Rot_lbc_eta_sig_r_eta_r2_22 = Rot_etsig_nu_r_21.*Rot_lcs_xi_nu_12 + Rot_etsig_nu_r_22.*Rot_lcs_xi_nu_22 + Rot_etsig_nu_r_23.*Rot_lcs_xi_nu_32;
Rot_lbc_eta_sig_r_eta_r2_23 = Rot_etsig_nu_r_21.*Rot_lcs_xi_nu_13 + Rot_etsig_nu_r_22.*Rot_lcs_xi_nu_23 + Rot_etsig_nu_r_23.*Rot_lcs_xi_nu_33;
%
Rot_lbc_eta_sig_r_sig_21 = Rot_etsig_nu_r_31.*Rot_lcs_xi_nu_11 + Rot_etsig_nu_r_32.*Rot_lcs_xi_nu_21 + Rot_etsig_nu_r_33.*Rot_lcs_xi_nu_31;
Rot_lbc_eta_sig_r_sig_22 = Rot_etsig_nu_r_31.*Rot_lcs_xi_nu_12 + Rot_etsig_nu_r_32.*Rot_lcs_xi_nu_22 + Rot_etsig_nu_r_33.*Rot_lcs_xi_nu_32;
Rot_lbc_eta_sig_r_sig_23 = Rot_etsig_nu_r_31.*Rot_lcs_xi_nu_13 + Rot_etsig_nu_r_32.*Rot_lcs_xi_nu_23 + Rot_etsig_nu_r_33.*Rot_lcs_xi_nu_33; 
%
sigma_r = [Rot_lbc_eta_sig_r_sig_21,Rot_lbc_eta_sig_r_sig_22,Rot_lbc_eta_sig_r_sig_23]; 
eta_new_r2 = [Rot_lbc_eta_sig_r_eta_r2_21,Rot_lbc_eta_sig_r_eta_r2_22,Rot_lbc_eta_sig_r_eta_r2_23];
eta_new_r1 = [Rot_lbc_eta_sig_r_eta_r1_11,Rot_lbc_eta_sig_r_eta_r1_12,Rot_lbc_eta_sig_r_eta_r1_13];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
J = [1 0;0 -1]; % Mirror symmetric matrixs
c_radius = 1; 
c_11 = 1/c_radius; %These are only true for sphere. % One needs to change this for the different geometry Ofcourse. 
c_12 = 0;
c_21 = 0;
c_22 = 1/c_radius;
C = [c_11,c_12,c_21,c_22]; % this is only true for a sphere, for a general object it needs to be calculate locally
% Calculating the inverse of the Theta_Matrix
Theta_eta_i = [Rot_etasig_nuxi_11, Rot_etasig_nuxi_12 ,Rot_etasig_nuxi_21, Rot_etasig_nuxi_22];
det_Theta_eta_i = Theta_eta_i(:,1).*Theta_eta_i(:,4) - Theta_eta_i(:,2).*Theta_eta_i(:,3);
inv_Theta_eta_i = [Theta_eta_i(:,4),-Theta_eta_i(:,2),-Theta_eta_i(:,3),Theta_eta_i(:,1)];
inv_Theta_eta_i = inv_Theta_eta_i./det_Theta_eta_i; % Inverse of Theta Matrix
%
inv_Theta_eta_i_11 = inv_Theta_eta_i(:,1);
inv_Theta_eta_i_12 = inv_Theta_eta_i(:,2);
inv_Theta_eta_i_21 = inv_Theta_eta_i(:,3);
inv_Theta_eta_i_22 = inv_Theta_eta_i(:,4);
%Calculating the Reflecetd Gaussian beam strucutre
Gamma_r_temp_11 = Gamma_new_i_Q_1 + 2*cos(theta_i_vect).*(inv_Theta_eta_i_11.^2*c_11 + inv_Theta_eta_i_11.*inv_Theta_eta_i_21*c_12 + inv_Theta_eta_i_11.*inv_Theta_eta_i_12 *c_21 + inv_Theta_eta_i_12.*inv_Theta_eta_i_21*c_22);
Gamma_r_temp_12 = -Gamma_new_i_Q_2 + 2*cos(theta_i_vect).*(-inv_Theta_eta_i_11.*inv_Theta_eta_i_12*c_11 - inv_Theta_eta_i_11.*inv_Theta_eta_i_22*c_12 - inv_Theta_eta_i_12.^2*c_21 - inv_Theta_eta_i_12.*inv_Theta_eta_i_22*c_22);
Gamma_r_temp_21 = -Gamma_new_i_Q_3 + 2*cos(theta_i_vect).*(-inv_Theta_eta_i_11.*inv_Theta_eta_i_21*c_11 - inv_Theta_eta_i_21.^2*c_12 - inv_Theta_eta_i_22.*inv_Theta_eta_i_11*c_21 - inv_Theta_eta_i_21.*inv_Theta_eta_i_22*c_22);
Gamma_r_temp_22 = Gamma_new_i_Q_4 + 2*cos(theta_i_vect).*(inv_Theta_eta_i_21.*inv_Theta_eta_i_12*c_11 + inv_Theta_eta_i_21.*inv_Theta_eta_i_22*c_12 + inv_Theta_eta_i_12.*inv_Theta_eta_i_22*c_21 + inv_Theta_eta_i_22.^2*c_22);
Gamma_r_temp = [Gamma_r_temp_11,Gamma_r_temp_12,Gamma_r_temp_21,Gamma_r_temp_22];
%
sig_r_Q = sigma_i - 2*(dot(sigma_i,nu_vect,2)).*nu_vect; %Computing the reflected field direction, the reflected field shall acts as the new incident GB for our algorithm
sig_r_Q = sig_r_Q./vecnorm(sig_r_Q,2,2);
%% Rotation of the transverse coordinate 
%Step 5: Rotation of the transverse coordinate 
theta_sig_r = acos(sigma_r(:,3));
phi_sig_r = mod(atan2(sigma_r(:,2),sigma_r(:,1)),2*pi);
eta_tilde_1r = [cos(theta_sig_r).*cos(phi_sig_r),cos(theta_sig_r).*sin(phi_sig_r),-sin(theta_sig_r)];
eta_tilde_2r = [-sin(phi_sig_r),cos(phi_sig_r), zeros(length(phi_sig_r),1)];
Phi_eta  = [dot(eta_tilde_1r,eta_new_r1),dot(eta_tilde_1r,eta_new_r2),dot(eta_tilde_2r,eta_new_r1),dot(eta_tilde_2r,eta_new_r2)];
%Gamma_tild_r = ((Phi_eta))*(Gamma_i)*(transpose(Phi_eta));
a_11 = dot(eta_tilde_1r,eta_new_r1,2); 
b_11 = dot(eta_tilde_1r,eta_new_r2,2);
c_11 = dot(eta_tilde_2r,eta_new_r1,2);
d_11 = dot(eta_tilde_2r,eta_new_r2,2);
Gamma_r_temp_1 = Gamma_r_temp(:,1); %Gamma_i_11
Gamma_r_temp_2 = Gamma_r_temp(:,2); %Gamma_i_12
Gamma_r_temp_3 = Gamma_r_temp(:,3); %Gamma_i_21
Gamma_r_temp_4 = Gamma_r_temp(:,4);%Gamma_i_22
% Roation of the Gamma_i_Q in the plane of incidence
Gamma_new_r_1 = a_11.^2.*Gamma_r_temp_1 + a_11.*b_11.*Gamma_r_temp_2 + a_11.*b_11.*Gamma_r_temp_3 + b_11.^2.*Gamma_r_temp_4;
Gamma_new_r_2 = a_11.*c_11.*Gamma_r_temp_1 + a_11.*d_11.*Gamma_r_temp_2 + c_11.*b_11.*Gamma_r_temp_3 + b_11.*d_11.*Gamma_r_temp_4;
Gamma_new_r_3 = a_11.*c_11.*Gamma_r_temp_1 + c_11.*b_11.*Gamma_r_temp_2 + a_11.*d_11.*Gamma_r_temp_3 + b_11.*d_11.*Gamma_r_temp_4;
Gamma_new_r_4 = c_11.^2.*Gamma_r_temp_1 + c_11.*d_11.*Gamma_r_temp_2 + c_11.*d_11.*Gamma_r_temp_3 + d_11.^2.*Gamma_r_temp_4;
Gamma_new_r = [Gamma_new_r_1,Gamma_new_r_2,Gamma_new_r_3,Gamma_new_r_4];% Rotation of the Gamma_i_Q in the plane of incidence
%% RCS Calculations in the far zone.
theta_obser = pi; %Observation angles in the far zone
phi_obser = pi; 
r_hat = [zeros(length(phi_sig_r),1) + sin(theta_obser)*cos(phi_obser),...
    zeros(length(phi_sig_r),1) + sin(theta_obser)*sin(phi_obser),...
    zeros(length(phi_sig_r),1) + cos(theta_obser)]; % Monostatic observation angle
theta = acos(r_hat(:,3));
r = 10000000000000*a;%Far Zone
phi = mod(atan2(r_hat(:,2),r_hat(:,1)),2*pi);
sig_dir = sin(theta_sig_r).*cos(phi_sig_r).*(r*sin(theta).*cos(phi)-Q_S(:,1)) + sin(theta_sig_r).*sin(phi_sig_r).*(r*sin(theta).*sin(phi)-Q_S(:,2))...
         + cos(theta_sig_r).*(r*cos(theta)-Q_S(:,3));
a_11 = sin(theta_obser).*cos(theta_sig_r).*cos(phi_obser-phi_sig_r)-cos(theta_obser).*sin(theta_sig_r);
b_11 = sin(theta_obser).*sin(phi_obser-phi_sig_r);
%
det_Gamma_new_r = Gamma_new_r(:,1).*Gamma_new_r(:,4) - Gamma_new_r(:,2).*Gamma_new_r(:,3);
inv_Gamma_new_r = [Gamma_new_r(:,4),-Gamma_new_r(:,2),-Gamma_new_r(:,3),Gamma_new_r(:,1)];
inv_Gamma_new_r = inv_Gamma_new_r./det_Gamma_new_r; 
%
inv_Gamma_new_r_11 = inv_Gamma_new_r(:,1);
inv_Gamma_new_r_12 = inv_Gamma_new_r(:,2);
inv_Gamma_new_r_21 = inv_Gamma_new_r(:,3);
inv_Gamma_new_r_22 = inv_Gamma_new_r(:,4);
% 
ff_GB_decay = exp(1j*k/2*(inv_Gamma_new_r_11.*a_11.^2+ inv_Gamma_new_r_12.*a_11.*b_11 + inv_Gamma_new_r_21.*a_11.*b_11 + inv_Gamma_new_r_22.*b_11.^2)); 
%Far field Gaussian Beam angular decay term
ff_phase_shift = exp(1j*k*dot(r_hat,Q_S,2)); % Far field fraunhofer phase shift term
Psi_thet_phi_mnst = A_i_Q.*sqrt(1./(det_Gamma_new_r)).*ff_GB_decay.*ff_phase_shift.*((sig_dir>0.85*r)==1);
Sph_mnst =  sum(Reconcoeff.*Psi_thet_phi_mnst.*exp(-1j*k*l_b_vect));
%Z_r(cnt1,14) = a_mu*Psi_thet_phi_mnst*exp(-1j*k*l_b);