%% Multiband analysis. 
function [k_subbands,RCS_array_multiband,time,N_arr_tot,beam_array_eps_tot,err_tot] = UWB_GBSM_Sph(k_min,k_max,del_k,a,b,z_explane,nu_max,theta_inc,phi_inc,theta_obser,phi_obser,swtch_scheme)
    b_bar = b;
    J = log2(k_max/k_min); % Total number of subbands in the entire frequency band.
    RCS_array_multiband = []; %Preallocate this variable, since it changes its structure in every loop 
    k_subbands = [];
    N_arr_tot = [];
    time = []; 
    beam_array_eps_tot = [];
    err_tot = [];
    %upsamp_sch = 0;% upsamp_sch = 0 for upsampling in x and 1 in upsampling in xi 
    upsamp_sch = mod(swtch_scheme,2); %Test
    for i = 1:1:J
        k_max = 2^i*k_min;
        k_temp_low = k_max/2;
       %k_iter = k_temp_low:2^i*del_k:2*k_temp_low;
        k_iter = linspace(k_temp_low,2*k_temp_low,2*del_k);
        if i == 1
            tic
            Z = Chapter4Expansionplane(a,b_bar,k_max,nu_max,z_explane,theta_inc,phi_inc);
            %[Z_next,RCS_array_oneoct] = beam_RCS_oneoct(Z,a,b_bar,k_max,nu_max,z_explane,theta_inc,phi_inc,theta_obser,phi_obser,del_k,upsamp_sch);
            [Z_next,RCS_array_oneoct,N_arr,beam_array_eps1_mnst_arr,err] = beam_RCS_oneoct(Z,a,b_bar,k_max,nu_max,z_explane,theta_inc,phi_inc,theta_obser,phi_obser,del_k,upsamp_sch);
            k_subbands = [k_subbands,k_iter];
            RCS_array_multiband = [RCS_array_multiband,RCS_array_oneoct];
            N_arr_tot = [N_arr_tot, N_arr];
            beam_array_eps_tot = [beam_array_eps_tot;beam_array_eps1_mnst_arr];
            err_tot = [err_tot; err];
            t = toc; 
            time = [time t];
        end
        if i > 1 %Beam expansion for the next lattice
            xi_bar = sqrt(2*pi*nu_max/(b_temp*k_max));
            Z_next(find((Z_next(:,1)*xi_bar).^2 + (Z_next(:,2)*xi_bar).^2 >1),:) = [];% Elimination of the beams propagating in evanescent mode
            theta_n = acos(sqrt(1-xi_bar^2*(Z_next(:,2).^2+Z_next(:,1).^2)));
            phi_n = mod(atan2(Z_next(:,2),Z_next(:,1)),2*pi);
            Z_next(:,5) = z_explane;
            Z_next(:,6) = theta_n;
            Z_next(:,7) = phi_n;
            Z_next(find((Z_next(:,1)*xi_bar).^2 + (Z_next(:,2)*xi_bar).^2 >1),:) = [];% Elimination of the beams propagating in evanescent mode
            theta_n = acos(sqrt(1-xi_bar^2*(Z_next(:,2).^2+Z_next(:,1).^2)));
            phi_n = mod(atan2(Z_next(:,2),Z_next(:,1)),2*pi);
            Z_next(:,5) = z_explane;
            Z_next(:,6) = theta_n;
            Z_next(:,7) = phi_n;
            tic
            %[Z_next,RCS_array_oneoct] = beam_RCS_oneoct(Z_next,a,b_bar,k_max,nu_max,z_explane,theta_inc,phi_inc,theta_obser,phi_obser,del_k,upsamp_sch);
            [Z_next,RCS_array_oneoct,N_arr,beam_array_eps1_mnst_arr,err] = beam_RCS_oneoct(Z_next,a,b_temp,k_max,nu_max,z_explane,theta_inc,phi_inc,theta_obser,phi_obser,del_k,upsamp_sch);
            k_subbands = [k_subbands,k_iter];
            RCS_array_multiband = [RCS_array_multiband,RCS_array_oneoct];
            N_arr_tot = [N_arr_tot, N_arr];
            beam_array_eps_tot = [beam_array_eps_tot;beam_array_eps1_mnst_arr];
            err_tot = [err_tot; err];
            t = toc;
            time = [time t];
        end
        %Beam expansion for the next lattice 
        switch swtch_scheme %Switch cases  to compute the RCS of the target by all 4 schemes successfully deployed 
         case 1
             b_temp = b_bar*(2^(-i));% X-sampling scheme
             upsamp_sch = 1; 
         case 2
             b_temp = b_bar*(2^(i));% X-sampling scheme
             upsamp_sch = 0; 
         case 3
             if mod(i+1,2) == 0
                 b_temp = b_bar/2;
                 upsamp_sch = 0;
             else
                 b_temp = b_bar;
                 upsamp_sch = 1;
             end
         case 4
             if mod(i+1,2) == 0
                 b_temp = 2*b_bar;
                 upsamp_sch = 1; 
             else
                 b_temp = b_bar;
                 upsamp_sch = 0; % x upsampling
             end
        end
    end
end


