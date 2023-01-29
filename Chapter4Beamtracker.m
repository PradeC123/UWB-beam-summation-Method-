%% Chapter 4 Part 3: Beam Tracking 
% This code is dedicated for track all the beams within the target domain
%% Input Parameters
% Z = Z(:,1:7) matrix that describes the frame parameter
% c_1 = Center of the Target 
% a = radius of the target
% nu_max = Overcompleteness Parameter 
%k_max = highest operating frequency of the band 
% k = Operating frequency 
% b = collimation distance of the ID-GB 
% z_explane = Location of the expansion plane
% a_mu = reconstruction coefficient at the frequency k
%% Note
% Please note that the whole function is frequency independent (k
% independent)
%% Code 
function [Beamtrack,Z,a_mu] = Chapter4Beamtracker(Z,c1,a,nu_max,k_max,b,z_explane)
    x_bar = sqrt(2*pi*nu_max*b/k_max);%Frequency independent lattice
    xi_bar = sqrt(2*pi*nu_max/(b*k_max));%Frequency Independent Lattice 
    Z_r = zeros(length(Z(:,1)),11);
    QS = zeros(length(Z(:,1)),3);
    theta_r = zeros(length(Z(:,1)),1);
    phi_r = zeros(length(Z(:,1)),1);
    l_b = zeros(length(Z(:,1)),1);
    hit = zeros(length(Z(:,1)),1);
    Gamma_r = zeros(length(Z(:,1)),4);
    A_Q = zeros(length(Z(:,1)),1);
    r1 = a;
    % %c2 = [-1*a,0,0];
    % %r2 = a; 
    % [X1,Y1,Z1] = sphere;
    % % r = 200;
    % % X2 = X1 * r;
    % % Y2 = Y1 * r;
    % % Z2 = Z1 * r;
    % % surf(X2,Y2,Z2);
    % hold on;
    for cnt1 = 1:1:length(Z(:,1))
        m_i1 = Z(cnt1,3);
        m_i2 = Z(cnt1,4); 
        X_0 = [m_i1*x_bar,m_i2*x_bar,z_explane];
        theta_i = Z(cnt1,6);
        phi_i = Z(cnt1,7);
        Gamma_i = [1/(1j*b*cos(theta_i)^2),0;0,1/(1j*b)];
        Gamma_i = (reshape(transpose(Gamma_i),[1,4]));
        %-----------------------------------------% 
        [QS(cnt1,:),theta_r(cnt1),phi_r(cnt1),l_b(cnt1),Gamma_r(cnt1,:),A_Q(cnt1),hit(cnt1)] = multi_intersection_3D(X_0,theta_i,phi_i,Gamma_i,c1,r1);%,c2,r2);   
    end
    Z_r(:,1:3) = QS(:,1:3);
    Z_r(:,4) = theta_r;
    Z_r(:,5) = phi_r;
    Z_r(:,6) = l_b;
    Z_r(:,7:10) = Gamma_r(:,1:4);
    Z_r(:,11) = A_Q; %FF beam profile
    Z_r(:,12) = hit;%(n,m) that hi ts the sphere
    Z_r(find(hit==0),:) = [];% Deleteing all the elements that corresponds to hit = 0
    Z(find(hit==0),:) = [];
    Beamtrack = Z_r;
    Z = Z;
    a_mu = Z(:,8);
end










