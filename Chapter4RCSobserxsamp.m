%% Chapter 4 Part 4: RCS-x-sampling 
% This code is dedicated to calculate the RCS of the target at any given
% observation point by an allsampling method for a UWB upto 4 octaves
%% Input Parameters
% Z = Z(:,1:7) matrix that describes the frame parameter
% c_1 = Center of the Target 
% theta_obser = observation point (Elevation angle)
% phi_obser = observation point (Azimuthal angle)
% a = Radius of the target
% a_mu = reconstruction coefficient at the frequency k
% k = Operating frequency
% eps1 = threshold of an epsilon of 10^-2
% eps2 = threshold of an epsilon of 10^-3
%% Code 
function [RCS_obser,beamseps1,beamseps2] = Chapter4RCSobserxsamp(Z,theta_obser,phi_obser,a,b,k_min,k_max,nu_max,z_b,c1,eps1,eps2)
    N = log2(k_max/k_min);%
    j = 1; %First band
    while(j<N)
        Z = Chapter4Expansionplane(a,b,2^(j+1)*k_min,nu_max,z_b);
        [a_mu,Z] = Chapter4Reconstructioncoeff(Z,theta_obser,phi_obser,nu_max,2^(j+1)*k_min,k,b,z_b);
        Chapter4Refbeammn(Z,2^(j+1)*k_min,k,0,-1,1,1,0.15);
        [Z_r,Z,a_mu] = Chapter4Beamtracker(Z,c1,a,nu_msx,2^(j+1)*k_min,b,z_b);
        [RCS_obser,beamseps1,beamseps2] = Chapter4RCSobser(Z_r,pi,pi,1,a_mu,k,eps1,eps2);
    end
  
end




