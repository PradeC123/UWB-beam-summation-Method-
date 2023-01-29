%% Validity of the GBSM region compared to the physical optics 
% This code's purpose is to obtain the k_min for which the physical optics
% results can be replaced by GBSM results. 
%% Output of the code 
% acc_arr : Tracks the error of the absolute value of the difference
% between the GO and the GBSM tolerance level is calculated
%% Code Input Parameters 
k_min = 2000;
k_max = 4000;
del_k = 50; %Jumps 
theta_obser = pi/3;
phi_obser = 0;
theta_inc = 0;
phi_inc = 0; 
nu_max = 0.15; 
a = 1;
b = 1; 
%% Code Functionality 
k_arr = k_min:del_k:k_max;
acc_arr = zeros(1,length(k_arr));
for i = 1:1:length(k_arr)
    k = k_arr(i);
    [scatt,RCS] = Errscatt(k,k_max,nu_max,a,b,theta_inc,phi_inc,theta_obser,phi_obser);
    u_PO = 1/2*exp(1j*2*k*1*cos(pi/2-theta_obser/2));
    accc_arr(i) = abs(u_PO-scatt)/abs(scatt);
    disp(i);
end 

