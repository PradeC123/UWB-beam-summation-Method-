% Date - 30/06/2020 (16:44) Israel Daylight Time

% Date - 30/06/2020 (20:11) Israel Daylight Time
% Updates: 1) Added a FF phase shift expression R_Q in the FF expression
%         
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

%  Psi_thet_phi_mnst is a 1X1 matrix, that describes the Monostatic FF
%  expression of the beam. (mnst is monostatic)
%-------------------------------------------------------------------------%
% Monostatic Conditions 

% \theti = \pi -\theta_i
% \phi = \pi + \phi_i % observation direction in the monostatic case
%-------------------------------------------------------------------------%
% atan2(y,x): atan2(Y,X) is the four quadrant arctangent of the elements of X and Y
%    such that -pi <= atan2(Y,X) <= pi. hence we need to use 
% alpha_eta .* (alpha_eta >= 0) + (alpha_eta + 2 * pi) .* (alpha_eta < 0),
%-------------------------------------------------------------------------%
function [Psi_thet_phi_mnst,ff_GB_decay,ff_phase_shift] = GB_FF_monost(Gamma_r,Rot_lbc_eta_sig_r,Amp_GB_r,theta_inc,phi_inc,k,R_Q)
% Reflected GB spherical angles
sigma_r = Rot_lbc_eta_sig_r(3,:); %reflected beam axis unit vector 
theta_sig_r = acos(sigma_r(3)); %reflected beam's elevation angle 
phi_sig_r = atan2(sigma_r(2),sigma_r(1)); %reflected beams's azimuthal angle
phi_sig_r = phi_sig_r .* (phi_sig_r >= 0) + (phi_sig_r + 2 * pi) .* (phi_sig_r < 0); 
% Beam transverse coordinate rotation
eta_tilde_1r = [cos(theta_sig_r)*cos(phi_sig_r),cos(theta_sig_r)*sin(phi_sig_r),-sin(theta_sig_r)];
eta_tilde_2r = [-sin(phi_sig_r),cos(phi_sig_r),0];
eta_1r = Rot_lbc_eta_sig_r(1,:);
eta_2r = Rot_lbc_eta_sig_r(2,:);
Phi_eta  = [ dot(eta_tilde_1r,eta_1r),dot(eta_tilde_1r,eta_2r);dot(eta_tilde_2r,eta_1r),dot(eta_tilde_2r,eta_2r)];
%n2_z = Rot_lbc_eta_sig_r(2,3);
%n1_z = Rot_lbc_eta_sig_r(1,3);
%alpha_eta = atan2(n2_z,n1_z);
%alpha_eta = alpha_eta .* (alpha_eta >= 0) + (alpha_eta + 2 * pi) .* (alpha_eta < 0);
%Phi_eta = [ cos(alpha_eta),sin(alpha_eta); -sin(alpha_eta),cos(alpha_eta)];
% Gamma_r rotation due to the beam coordinates rotation
Gamma_tild_r = ((Phi_eta))*Gamma_r*(transpose(Phi_eta));
%FF(Far field) expression, the monostatic case.
theta =  pi - theta_inc;
phi = pi + phi_inc;
omega = [theta-theta_sig_r; sin(theta)*(phi-phi_sig_r)];% (Far field angular coordinate system)this applies for the monostatic conditions only
%omega = [sin(theta)*cos(phi)-cos(phi_sig_r)*sin(theta_sig_r);sin(theta)*sin(phi)-sin(phi_sig_r)*sin(theta_sig_r)];
%omega = [sin(theta_sig_r);0];
r_hat = [sin(theta)*cos(phi),sin(theta)*sin(phi),cos(theta)]; % Monostatic observation angle
ff_GB_decay = exp(1j*k/2*transpose(omega)*(inv(Gamma_tild_r))*omega);% Far field Gaussian Beam angular decay term
ff_phase_shift = exp(1j*k*dot(r_hat,R_Q)); % Far field fraunhofer phase shift term
Psi_thet_phi_mnst = Amp_GB_r/(sqrt(det(Gamma_tild_r)))*ff_GB_decay*ff_phase_shift;
end 