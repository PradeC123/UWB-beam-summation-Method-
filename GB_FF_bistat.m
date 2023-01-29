% Date - 30/06/2020 (19:10) Israel Daylight Time
%-------------------------------------------------------------------------%
% Please read the document "Reflection of a General Gaussian Beam by a 3D Curved Surface" 
% to undertand the working of the function. Specifically "2. Far field
% expressions for a general Gaussian Beam"
%-------------------------------------------------------------------------%
% Rot_lbc_eta_sig_r = Rotation of the local beam coordinates
% (\sigma_r,\4\eta_r) describes the projection of the reflected beam
% coordinates on the global coordinate systems. 
%-------------------------------------------------------------------------%
% INPUT of the function
% Gamma_r is the 2X2 complex symmetrical matrix of the the reflected beam

% Rot_lbc_eta_sig_r is the 3X3 projection matrix that describes the projection
% of the reflected beam coordinates on the global coordinate systems. 

% Amp_GB_r is the 1X1 matrix, that describes the amplitude of the reflected
% beam.

% theta_inc is the 1X1 matrix, the describes the initial elevation angle of
% the incidence source GB beam

% phi_inc is the 1X1 matrix, that describes the inital azimuthal angle of
% incidence source GB beam 

% k is the 1X1 matrix, that describes the wavenumber of the reflected beam

% R_Q is the 1X3 matrix, that describes the last hit point to the curve,
% before going to the far field. which will be accounted in the FF beam expression as the the phase
% shift in the far zone region
%-------------------------------------------------------------------------%
% OUTPUT of the function 

%  Psi_thet_phi_bst is a 1X1 matrix, that describes the Monostatic FF
%  expression of the beam. (bst is bistatic)
%-------------------------------------------------------------------------%
% atan2(y,x): atan2(Y,X) is the four quadrant arctangent of the elements of X and Y
%    such that -pi <= atan2(Y,X) <= pi. hence we need to use 
% alpha_eta .* (alpha_eta >= 0) + (alpha_eta + 2 * pi) .* (alpha_eta < 0),
%-------------------------------------------------------------------------%
function [Psi_thet_phi_bstat] = GB_FF_bistat(Gamma_r,Rot_lbc_eta_sig_r,Amp_GB_r,k,R_Q)
% Meshgrid (Theta,Phi)
theta_glb = 0:0.001:pi; %glb is global.
phi_glb = 0:0.001:2*pi;
[Theta_glb,Phi_glb] = meshgrid(theta_glb,phi_glb);
Psi_thet_phi_bstat = zeros(length(phi_glb),length(theta_glb));
% Reflected Beam spherical angles
sigma_r = Rot_lbc_eta_sig_r(3,:); %reflected beam axis unit vector 
theta_sig_r = acos(sigma_r(3)); %reflected beam's elevation angle 
phi_sig_r = atan2(sigma_r(2),sigma_r(1)); %reflected beams's azimuthal angle
phi_sig_r = phi_sig_r .* (phi_sig_r >= 0) + (phi_sig_r + 2 * pi) .* (phi_sig_r < 0); 
% Beam transverse coordinate rotation
n2_z = Rot_lbc_eta_sig_r(2,3);
n1_z = Rot_lbc_eta_sig_r(1,3);
alpha_eta = atan2(n2_z,n1_z);
alpha_eta = alpha_eta .* (alpha_eta >= 0) + (alpha_eta + 2 * pi) .* (alpha_eta < 0);
Phi_eta = [ cos(alpha_eta),sin(alpha_eta); -sin(alpha_eta),cos(alpha_eta)];
% Gamma_r rotation due to the beam coordinates rotation
Gamma_tild_r = transpose(inv(Phi_eta))*Gamma_r*(inv(Phi_eta));
%FF expression
r_hat_dt_R_Q = sin(Theta_glb).*cos(Phi_glb)*R_Q(1)+ sin(Theta_glb).*sin(Phi_glb)*R_Q(2) + cos(Theta_glb)*R_Q(3);
ff_phase = exp(1j*k*r_hat_dt_R_Q);
Gamma_ff = imag(inv(Gamma_tild_r));
Gamma_ff_11 = Gamma_ff(1,1);
Gamma_ff_12 = Gamma_ff(1,2);
Gamma_ff_21 = Gamma_ff(2,1);
Gamma_ff_22 = Gamma_ff(2,2);
ff_beam_amp_quad_fm = Gamma_ff_11*(Theta_glb-theta_sig_r).^2 + ...
    Gamma_ff_12*sin(theta_sig_r)*(Theta_glb-theta_sig_r).*(Phi_glb - phi_sig_r) + ...
    Gamma_ff_21*sin(theta_sig_r)*(Theta_glb-theta_sig_r).*(Phi_glb - phi_sig_r) + ...
    Gamma_ff_22*sin(theta_sig_r).^2*(Phi_glb - phi_sig_r).^2 ;
Psi_thet_phi_bstat = Amp_GB_r/(sqrt(det(Gamma_tild_r))) * exp(-k/2*ff_beam_amp_quad_fm).*ff_phase;
end 