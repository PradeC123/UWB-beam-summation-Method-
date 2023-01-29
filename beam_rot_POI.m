% Date - 29/06/2020 (21:47) Israel Daylight Time
%-------------------------------------------------------------------------%
% Please read the document "Computation of the Local GB coordinates in the Plane of Incidence" 
% to undertand the working of the function
%-------------------------------------------------------------------------%
% Rot_lbc_eta_sig_new = Rotation of the local beam coordinates (\sigma,\4eta)
% beam_rot_POI = beam rotation in the Plane of Incidence

% Rot_lcs_xi_nu = Rotation of the local coordinate system in \4\xi and \nu
%-------------------------------------------------------------------------%
% INPUT of the function
% Gamma_i_int is the 2X2 complex symmetrical matrix of the incidence beam
% at the point of emergence of the beam.

% Rot_etasig_xyz is the 3X3 projection matrix that describes the projection of..
% the local beam coordinates on the global coordinate systems.

% nu_vect is the 1X3 Matrix that descrbibes the unit vector in the  direction of the beam propagation
%-------------------------------------------------------------------------%
% OUTPUT of the function 
% Gamma_new_i_int is the 2X2 complex symmetrical matrix of the the new
% rotated incidence beam at the point of the emergence of the beam.

% Rot_lbc_eta_sig_new is the 3X3 projection matrix that describes the projection
% of the new rotated local beam coordinates on the global coordinate systems.

% Rot_lcs_xi_nu is the 3X3 projection matrix that describes the projection
% of the local surface coordinates (\xi_1,\xi_2,\nu) on the global
% coordinate systems.(This coordinate system will be useful for the
% computation of the reflected beam (Please look GB_Refl.m). 
%-------------------------------------------------------------------------
%-------------------------------------------------------------------------%
function [Gamma_new_Q_inc,Rot_etasig_new_inc,Rot_Q_xinu_xyz] = beam_rot_POI(Gamma_i_init, Rot_etasig_xyz,nu_vect,l_b)
% sigma_i = Rot_etasig_xyz(:,7:9);
sigma_i = Rot_etasig_xyz(3,:);
xi_1 = cross(nu_vect,-sigma_i,2);
if vecnorm(xi_1,2,2)==0 %This case implies Plane of incidence doesnot exists and nu_nect and sigma are colinear.
    xi_1 = Rot_etasig_xyz(1,:);
    xi_2 = cross(nu_vect,xi_1,2);
else
    xi_1 = xi_1./vecnorm(xi_1,2,2);%normalization of the vector
    xi_2 = cross(nu_vect,xi_1,2);
    xi_2 = xi_2./vecnorm(xi_2,2,2);
end
% Choosing the desired beam coordinate system such that, eta_new_1 is in
% the perpendicular direction of the plane of incidence and eta_new_ is in the plane of
% incendence
eta_new_1 = xi_1;
eta_new_2 = cross(sigma_i,eta_new_1,2);
eta_new_2 = eta_new_2./vecnorm(eta_new_2,2,2);
eta_i_1 = Rot_etasig_xyz(1,:);
eta_i_2 = Rot_etasig_xyz(2,:);
% The matrix that rotates the plane of the transverse coordinate system to
% the new desired plane such that it satisfies the specified conditions as
% mentioned above.
theta_eta_12 = [ dot(eta_new_1,eta_i_1,2), dot(eta_new_1,eta_i_2,2); ...
    dot(eta_new_2,eta_i_1,2), dot(eta_new_2,eta_i_2,2)];
Gamma_i_Q = inv(inv(Gamma_i_init)+l_b*eye(2));
Gamma_new_Q_inc = (theta_eta_12)*Gamma_i_Q*(transpose(theta_eta_12)); 
Rot_etasig_new_inc = [eta_new_1;eta_new_2;sigma_i];
Rot_Q_xinu_xyz = [xi_1;xi_2;nu_vect];
%Rot_etasig_new_inc = cat(2,eta_new_1,eta_new_2,sigma_i);
%Rot_Q_xinu_xyz = cat(2,xi_1,xi_2,nu_vect);
end 
