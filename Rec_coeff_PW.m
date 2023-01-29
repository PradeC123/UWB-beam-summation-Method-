% Date - 02/09/2020 (18:42) Israel Daylight Time (Started on)
% 05/09/2020(19:56,Saturday) 
% This function is only meant for calculating the reconstuction coefficient
% for a target incident my plane waves.
function [a_mu] = Rec_coeff_PW(theta_0,phi_0,nu,b,k,m1,m2,n1,n2,z_b)
 xi_bar = sqrt(2*pi*nu/(k*b)); %spectral translation 
 x_bar = sqrt(2*pi*b*nu/k); %Spatial Translation
 a_mu = nu^2*2*exp(-k*b*(sin(theta_0)*cos(phi_0)-n1*xi_bar).^2/2).*...
 exp(-k*b*(sin(theta_0)*sin(phi_0)-n2*xi_bar).^2/2).*...
 exp(-1j*k*sin(theta_0)*cos(phi_0).*m1*x_bar).*exp(-1j*k*sin(theta_0)*sin(phi_0).*m2*x_bar).*exp(-1j*k*z_b*cos(theta_0));
end

