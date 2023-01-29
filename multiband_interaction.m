%% Intersection of a beam lattice with a Sphere (Frequency Independent)
function [Q_S,nu_vect,l_b_vect,Reconcoeff,Z]  = targetbeam_interaction(Z,b,nu_max,k_max,z_explane,c_1,r_1,Reconcoeff)
theta_n = Z(:,6);
phi_n = Z(:,7);
x_bar = sqrt(2*pi*nu_max*b/k_max);
xi_bar = sqrt(2*pi*nu_max/(b*k_max));
sigma_i = [sin(theta_n).*cos(phi_n),sin(theta_n).*sin(phi_n),cos(theta_n)];
X_0 = [Z(:,3)*x_bar,Z(:,4)*x_bar, zeros(length(Z(:,3)),1)+z_explane];
Delta_1 = dot(sigma_i,X_0-c_1,2).^2 - vecnorm(X_0-c_1,2,2).^2 + r_1^2;
Z(find(Delta_1<0),:) = []; %Removing the beams that doesnot hit the sphere
Reconcoeff(find(Delta_1<0),:) = [];
Delta_1(find(Delta_1<0),:) = []; %Removing the beams that doesnot hit the sphere
X_0(find(Delta_1<0),:) = []; %Removing the beams that doesnot hit the sphere
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
end

