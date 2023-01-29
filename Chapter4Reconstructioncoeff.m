%% Chapter 4 Part 2: Reconstruction Coeff
% This code is dedicated for calculating the reconstruction coefficient 
%% Input Parameters
% Z = Z(:,1:7) matrix that describes the frame parameter (Must be
% calculated via Chapter 4 Part 1: Expansion plane. 
% k = operating frequency 
% theta_inc = Initial angle of incidence (Elevation)
% phi_inc = Initial angle of incidence (Azimuthal) 
% k_max = highest frequency in the frequency band 
% k = frequency of the operating band 
% b = collimation distance of the ID-GB
% z_explane = location of the incident expansion  plane.
% nu_max = Overcompleteness parameter
% Note: Analysis only valid if the initial wave is plane wave 
%% Code 
function [Reconcoeff,Z] = Chapter4Reconstructioncoeff(Z,theta_inc,phi_inc,nu_max,k_max,k,b,z_explane)
    x_bar = sqrt(2*pi*nu_max*b/k_max);%Frequency independent lattice
    xi_bar = sqrt(2*pi*nu_max/(b*k_max));
    freq_indep = b*(sin(theta_inc)*cos(phi_inc)-Z(:,1)*xi_bar).^2/2 + ...
        b*(sin(theta_inc)*sin(phi_inc)-Z(:,2)*xi_bar).^2/2 + ...
        + 1j*sin(theta_inc)*cos(phi_inc)*Z(:,3)*x_bar + 1j*sin(theta_inc)*sin(phi_inc)*Z(:,4)*x_bar + ...
        + 1j*z_explane*cos(theta_inc); %Frequency independent terms in a_mu
    Reconcoeff = nu_max^2*2*(k/k_max).^2.*exp(-k.*freq_indep);  %Chapter 2, equation numbrt (2.25a) and (2.25b) 
    %Reconcoeff = nu_max^2*2*(k/k_max)^2*exp(-k*b*(sin(theta_inc)*cos(phi_inc)-Z(:,1)*xi_bar).^2/2).*exp(-k*b*(sin(theta_inc)*sin(phi_inc)-Z(:,2)*xi_bar).^2/2).*...
    %exp(-1j*k*sin(theta_inc)*cos(phi_inc)*Z(:,3)*x_bar).*exp(-1j*k*sin(theta_inc)*sin(phi_inc)*Z(:,4)*x_bar).*exp(-1j*k*z_explane*cos(theta_inc)); 
    if length(k) == 1 %handling the mutliband case
        Z(:,8) = Reconcoeff;
        Z(abs(Reconcoeff)<10^-8,:)=[];
        Reconcoeff = Z(:,8); 
    end
end