% This code is dedicated to calculate the RCS of the target at any given
% observation point for all the schemes provided4
%% Input Parameters
% Z = Z(:,1:7) matrix that describes the frame parameter
% c_1 = Center of the Target 
% theta_obser = observation point (Elevation angle)
% phi_obser = observation point (Azimuthal angle)
% a = Radius of the target
% a_mu = reconstruction coefficient at the frequency k
% k_min = Smallest Operating frequency of the band 
% k_max = Largest operating frequency of the band
% nu_bar = Largest overcompleteness parameter
% b_bar = Collimation distance according to the scheme  
% z_b = Location of the plane of expansion
% switch cases to activate the scheme 
% 1) X sampling scheme 2) Xi sampling scheme 3) Interlaced XXi sampling scheme 4)
% Interlaced XiX sampling scheme 
% eps1 = thresholding level 1 
% eps2 = thresho;ding level 2
%%%%%%%%%%%%%%%%%%%
%% Code 
function [RCS_array,beamstrack_array] = Chapter4RCSUWB(k_min,k_max,del_k,theta_inc,phi_inc,theta_obs,phi_obs,z_b,c_1,a,b_bar,nu_bar,eps1,eps2,swtch_scheme)
    J = log2(k_max/k_min); % Total number of subbands in the entire frequency band.
    % Code of x-sampling 
    k = k_min:del_k:k_max; %Entire frequency band
    k_subbands = []; % Storing all the octave-subbands in the whole frequency band 
    k_rlcnt = 0; % Global count of the total number of elements in frequencies in the band. 
    RCS_array = zeros(1,length(k));
    beamstrack_array = zeros(1,length(k));
    k_temp_low = k_min; % minimum frequecny of the band
    b_temp = b_bar; % temp value of b, since it changes with the scheme
    for i = 1:1:J
         Z = Chapter4Expansionplane(a,b_temp,2*k_temp_low,nu_bar,z_b);
         [a_mu,Z] = Chapter4Reconstructioncoeff(Z,theta_inc,phi_inc,nu_bar,2*k_temp_low,2*k_temp_low,b_temp,z_b);
         [Z_r,Z,a_mu] = Chapter4Beamtracker(Z,c_1,a,nu_bar,2*k_temp_low,b_temp,z_b);
         k_iter = k_temp_low:del_k:2*k_temp_low; % Single octave
         k_subbands = [k_subbands,k_iter];
         for k_cnt = 1:1:length(k_iter)% Tracing all the beams at the highest frequency of the octave
              k_rlcnt = k_rlcnt + 1;
              [a_mu,Z] = Chapter4Reconstructioncoeff(Z,theta_inc,phi_inc,nu_bar,2*k_temp_low,k_iter(k_cnt),b_temp,z_b);
              [RCS_obser,beamseps1,beamseps2] = Chapter4RCSobser(Z_r,theta_obs,phi_obs,a,a_mu,k_iter(k_cnt),eps1,eps2);
              RCS_array(k_rlcnt) = RCS_obser;
              beamstrack_array(k_rlcnt) = beamseps1;
         end 
         switch swtch_scheme %Switch cases  to compute the RCS of the target by all 4 schemes successfully deployed 
             case 1
                 b_temp = b_bar*(2^(-i));% X-sampling scheme
             case 2
                 b_temp = b_bar*(2^(i));% X-sampling scheme
             case 3
                 if mod(i+1,2) == 0
                     b_temp = b_temp/2;
                 else
                     b_temp = b_temp;
                 end
             case 4
                 if mod(i+1,2) == 0
                     b_temp = 2*b_temp;
                 else
                     b_temp = b_temp;
                 end
         end 
         k_temp_low = 2*k_temp_low;% Transition towards the next sub-octave in the frequency band.  
    end
end
