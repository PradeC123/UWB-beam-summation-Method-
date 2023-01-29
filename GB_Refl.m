% Date - 29/06/2020 (19:51) Israel Daylight Time

% Latest updates: 1) added an amplitude parameter in the reflected beam 
%                 2) Added the reference equation number to be refered from
%                 the document specified below.
% Date -30/06/2020 (15:44) Israel Daylight Time
% 
%-------------------------------------------------------------------------%
% Latest updates: 1) added the phase parameter in the amplitude term
                     %an amplitude parameter in the reflected beam
                     %{exp(-1j*k*l_b_GB)}
% Date - 12/11/2020 (03:50) Israel Daylight Time

%-------------------------------------------------------------------------%
% Please read the document "Reflection of a General Gaussian Beam by a 3D Curved Surface" 
% to undertand the working of the function, Specifically the reflection
% part EQ no. (4),(5),(6),(7),(10),(17).
%-------------------------------------------------------------------------%
% Rot_lbc_eta_sig = Rotation of the local beam coordinates (\sigma,\4\eta)
% beam_rot_POI = beam rotation in the Plane of Incidence
% Rot_lcs_xi_nu = Rotation, local coordinate system (\nu,\4\xi)
%-------------------------------------------------------------------------%
% INPUT of the function

% Gamma_new_i_int is the 2X2 complex symmetrical matrix of the the new
% rotated incidence beam at the point of its emergence

% Rot_lbc_eta_sig_new_int is the 3X3 projection matrix that describes the projection
% of the new rotated local beam coordinates on the global coordinate systems.

% Rot_lcs_xi_nu is the 3X3 projection matrix that describes the projection
% of the local surface coordinates (\xi_1,\xi_2,\nu) on the global
% coordinate systems.(This coordinate system will be useful for the
% computation of the reflected beam.

% l_b is the distance propagated by the beam from the emergence point to
% the point where the beam hits the interface Q, (in our case a sphere). 

% r is the 1X1 matrix, which is the radius of the sphere.

% k is the wavenumber of the beam
%-------------------------------------------------------------------------%
% OUTPUT of the function 
% Gamma_r is the 2X2 complex symmetrical matrix of the the reflected beam

% Rot_lbc_eta_sig_r is the 3X3 projection matrix that describes the projection
% of the reflected beam coordinates on the global coordinate systems. 

% Amp_GB_r is the 1X1 matrix, that describes the amplitude of the reflected
% beam.
%-------------------------------------------------------------------------%
function [Gamma_r,Rot_lbc_eta_sig_r,Amp_GB_r] = GB_Refl(Gamma_i_init,Gamma_new_i_int,l_b_GB,Rot_lbc_eta_sig_new_int,Rot_lcs_xi_nu,r)
% eta_i_1, eta_i_2 and sigma_i, vectors obtained by rotating the old
% coordinate system into the new coordinate system, such that eta_i_1 is
% perpendicular to the plane of incidence and eta_i_2 is in the plane of
% incidence.
eta_new_i_1 = Rot_lbc_eta_sig_new_int(1,:);
eta_new_i_2 = Rot_lbc_eta_sig_new_int(2,:);
sigma_i = Rot_lbc_eta_sig_new_int(3,:);
% xi_1, xi_2 and nu are the local coordinate on the sphere, at the point Q
xi_1 = Rot_lcs_xi_nu(1,:);
xi_2 = Rot_lcs_xi_nu(2,:);
nu = Rot_lcs_xi_nu(3,:);
theta_i = acos(dot(-sigma_i,nu));% the angle between the sigma_i and nu, the initial angle of incidence.
% Matrix relation between (eta_i_1, eta_i_2, sigma_i) and (xi_1, xi_2, nu)
Rot_etsig_nu_i = [dot(eta_new_i_1,xi_1,2), dot(eta_new_i_1,xi_2,2), dot(eta_new_i_1,nu,2);...
    dot(eta_new_i_2,xi_1,2),dot(eta_new_i_2,xi_2,2),dot(eta_new_i_2,nu,2);...
    dot(sigma_i,xi_1,2),dot(sigma_i,xi_2,2), dot(sigma_i,nu,2)];
% Projection matrix of (eta_i_1,eta_i_2) on (xi_1,xi_2)
 Theta_eta_i = Rot_etsig_nu_i([1,2],[1,2]);
% Matrix relation between (eta_r_1, eta_r_2, sigma_i) and (xi_1, xi_2, nu) 
Rot_etsig_nu_r = [dot(eta_new_i_1,xi_1,2), dot(eta_new_i_1,xi_2,2), dot(eta_new_i_1,nu,2);...
    dot(eta_new_i_2,xi_1,2) -dot(eta_new_i_2,xi_2,2),dot(eta_new_i_2,nu,2);...
    dot(sigma_i,xi_1,2),dot(sigma_i,xi_2,2), -dot(sigma_i,nu,2)];
%Reflected beam coordinates with respect to the glbl cordinates.
Rot_lbc_eta_sig_r = Rot_etsig_nu_r * Rot_lcs_xi_nu; 
% The behaviour of the Gamma(sigma_i) along the sigma_i axis
% Analysis for the Gamma_r for the reflected beam
J = [1 0;0 -1]; % Mirror symmetric matrix s
C = [1/r,0;0,1/r]; % this is only true for a sphere, for a general object it needs to be calculate locally.s
Gamma_r = J*(Gamma_new_i_int)*J + 2*cos(theta_i)*J*transpose(inv(Theta_eta_i))*C*(inv(Theta_eta_i))*J;
Amp_GB_r = sqrt(det(Gamma_new_i_int)/det(Gamma_i_init));
end


