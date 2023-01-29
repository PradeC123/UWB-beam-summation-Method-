%% Chapter 4 Part 1: Expansion Plane
% This code is dedicated for location of the expansion for UWB and
% establishing an expansion plane that shoots ID-GB's from this plane 
%% Input Parameters
% theta_inc = Initial angle of incidence(Elevation)
% phi_inc = Initial angle of incidence(Azimuthal) 
% a = Radius of the target 
% b = Collimatiion distance of the ID-GB
% k_max = Highest frequency in the frequency band
% nu_max = Overcompleteness Parameter 
% z_explane = Location of the expansion plane
%% Code 
function Expplane = Chapter4Expansionplane(a,b,k_max,nu_max,z_explane,theta_inc,phi_inc)
    x_bar = sqrt(2*pi*nu_max*b/k_max);%Frequency independent lattice
    xi_bar = sqrt(2*pi*nu_max/(b*k_max)); % Chapter 2 equation (2.14) in the thesis
    z_b = z_explane; %location of the expansion plane  
    %% Lattice expansion indices.
    %L = x1_max; %the width of the square aperture in x1 axis. Onl relevant for truncation of the beams 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Respective indices in m_1,m_2,n_1,n_2
    M1 = floor(2*a/x_bar); %Parameters M1, M2 should be chooses such that the expansion plane spans the whole target domain  
    M2 = floor(2*a/x_bar);
    m_1 = -M1:1:M1;
    m_2 = -M2:1:M2;
    n_1max = round(1/xi_bar*sin(theta_inc)*cos(phi_inc));
    n_2max = round(1/xi_bar*sin(theta_inc)*sin(phi_inc));
    C = 144; %Optimization parameter
    del_n = round(sqrt(C/(pi*nu_max)));
    n_1 = n_1max - del_n:1:n_1max + del_n;
    n_2 = n_2max -del_n:1:n_2max + del_n;
    N1 = floor(sin(45*pi/180)/(xi_bar));
    N2 = floor(sin(45*pi/180)/(xi_bar));
    %N1 = 0;
    %N2 = 0;
    %n_1 = -N1:1:N1;
    %n_2 = -N2:1:N2;
    %n_1 = 1;
    %n_2 = 1;
    [N_1,N_2,M_1,M_2] = ndgrid(n_1,n_2,m_1,m_2);
    Z = [N_1(:),N_2(:),M_1(:),M_2(:)];% Generates all the possible combination of n1,n2,m1,m2
    Z((Z(:,1)*xi_bar).^2 + (Z(:,2)*xi_bar).^2 >1,:) = [];% Elimination of the beams propagating in evanescent mode
    %Z(find((Z(:,1)*xi_bar).^2 + (Z(:,2)*xi_bar).^2 >1),:) = [];% Elimination of the beams propagating in evanescent mode
    theta_n = acos(sqrt(1-xi_bar^2*(Z(:,2).^2+Z(:,1).^2)));
    phi_n = mod(atan2(Z(:,2),Z(:,1)),2*pi);
    Z(:,5) = z_b; 
    Z(:,6) = theta_n;
    Z(:,7) = phi_n;
    Expplane = Z;
end
