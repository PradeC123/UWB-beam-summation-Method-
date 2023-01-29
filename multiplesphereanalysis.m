%% Date: 10/09/2021 Israel Daytime 12:20
%%PO Results for monostatic RCS \theta = pi for Multisphere analysis
%% Initial Input Parameters 
% a is the radius of the sphere
% a+d is the distance from the origin to the sphere center
% theta_i is the initial angle of incidence
% rho_i_1 is the radius of the curvature in the plane of incidence
% rho_i_2 is the radius of the curvature in the perp. to the POI
%% Code Input Parameters 
k_arr_min = 1000;
k_arr_max = 1001;
del_k = 0.00001;
a = 1000;
theta_i = 0;
d = 1000;
%% Code Functionality
rho_i_1 = (a*cos(theta_i))/2;
rho_i_2 = a/(2*cos(theta_i)); 
k_arr = k_arr_min:del_k:k_arr_max;
u_tot_1  = 2*sqrt(rho_i_1*rho_i_2).*exp(1j*2*k_arr*a); %Single reflection backscattering 
s_0 = 2*(a-a/sqrt(2)+d);
rho_i_11 = a/(2*sqrt(2)); rho_i_22 = a/sqrt(2);
u_tot_2 = 2*sqrt(rho_i_1*rho_i_2).*exp(1j*sqrt(2)*k_arr*a).*exp(-1j*k_arr*s_0)/(sqrt(1+s_0/rho_i_11)*sqrt(1+s_0/rho_i_22));
% GO interpretation for multiple reflection back scattering
%% Plotting the results
plot(k_arr,4*pi*abs(u_tot_1+u_tot_2).^2);