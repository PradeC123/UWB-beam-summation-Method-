function [Z_next,RCS_array_oneoct,N_arr,beam_array_eps1_mnst_arr,err] = beam_RCS_oneoct(Z,a,b,k_max,nu_max,z_explane,theta_inc,phi_inc,theta_obser,phi_obser,del_k,upsamp_sch)
    %% Parameters of the sphere
    c_1 = [0,0,0]; %The location of the center of the sphere 
    r_1 = a; % The radius of the sphere 
    %Tar_surf = [c_1,r_1];%Target surface
    del = 2*1e9;
    mini_batch = floor(length(Z)/del);
    sum_sph_mnst = 0; 
    Z_i_mini_batches = [];
    Z_r_mini_batches = [];
    N_arr = [];
    beam_array_eps1_mnst_arr = [];
    err = [];
    k = k_max; % Calculating the Rcs for the hightest frequency of the band
    if mini_batch == 0 
        Z_i_mini = Z;
        [Reconcoeff,Z_i_mini] = Chapter4Reconstructioncoeff(Z_i_mini,theta_inc,phi_inc,nu_max,k_max,k,b,z_explane);
        [Q_S,nu_vect,l_b_vect,Reconcoeff,Z_i_mini]  = targetbeam_interaction(Z_i_mini,b,nu_max,k_max,z_explane,c_1,r_1,Reconcoeff);
        [Gamma_new_i_Q_Mat,A_i_Q,xi_1,xi_2,eta_new_i1,eta_new_i2,sigma_i] =  Gamma_i_rot(Z_i_mini,l_b_vect,b,k_max,nu_max,nu_vect);
        [Gamma_r_mat,eta_new_r1,eta_new_r2,sigma_r] =  GBReflection(xi_1,xi_2,nu_vect,eta_new_i1,eta_new_i2,sigma_i,Gamma_new_i_Q_Mat,l_b_vect);
        [Gamma_new_r_mat,Gamma_new_r,theta_sig_r,phi_sig_r] = N2FF_Gammar_rot(Gamma_r_mat,eta_new_r1,eta_new_r2,sigma_r);
        Z_r = [Q_S,theta_sig_r,phi_sig_r,l_b_vect,A_i_Q,Gamma_new_r]; 
        [sum_sph_mnst, Z_r] =  N2FF_GBSM(Gamma_r_mat,Gamma_new_r,l_b_vect,theta_obser,phi_obser,phi_sig_r,theta_sig_r,Q_S,k,A_i_Q,Reconcoeff,Z_r);
        Z_r_mini_batches = [Z_r_mini_batches;Z_r]; 
        Z_i_mini_batches = [Z_i_mini_batches;Z_i_mini];
        Z_beams = [Z_i_mini,Z_r(:,12)];    
        disp(4*abs(sum_sph_mnst).^2)
    else 
        for i_cnt = 1:1:mini_batch
            Z_i_mini = Z((i_cnt-1)*del+1:i_cnt*del,:);
            [Reconcoeff,Z_i_mini] = Chapter4Reconstructioncoeff(Z_i_mini,theta_inc,phi_inc,nu_max,k_max,k,b,z_explane);
            [Q_S,nu_vect,l_b_vect,Reconcoeff,Z_i_mini]  = targetbeam_interaction(Z_i_mini,b,nu_max,k_max,z_explane,c_1,r_1,Reconcoeff);
            [Gamma_new_i_Q_Mat,A_i_Q,xi_1,xi_2,eta_new_i1,eta_new_i2,sigma_i] =  Gamma_i_rot(Z_i_mini,l_b_vect,b,k_max,nu_max,nu_vect);
            [Gamma_r_mat,eta_new_r1,eta_new_r2,sigma_r] =  GBReflection(xi_1,xi_2,nu_vect,eta_new_i1,eta_new_i2,sigma_i,Gamma_new_i_Q_Mat,l_b_vect);
            [Gamma_new_r_mat,Gamma_new_r,theta_sig_r,phi_sig_r] = N2FF_Gammar_rot(Gamma_r_mat,eta_new_r1,eta_new_r2,sigma_r);
            Z_r = [Q_S,theta_sig_r,phi_sig_r,l_b_vect,A_i_Q,Gamma_new_r];
            [sum_sph_mnst_mini, Z_r_mini] = N2FF_GBSM(Gamma_new_r_mat,Gamma_new_r,l_b_vect,theta_obser,phi_obser,phi_sig_r,theta_sig_r,Q_S,k,A_i_Q,Reconcoeff,Z_r);
            Z_r_mini_batches = [Z_r_mini_batches;Z_r_mini];
            Z_i_mini_batches = [Z_i_mini_batches;Z_i_mini]; 
            sum_sph_mnst = sum_sph_mnst + sum_sph_mnst_mini;
            disp(i_cnt)
        end
        Z_i_mini = Z(i_cnt*del+1:length(Z),:);
        [Reconcoeff,Z_i_mini] = Chapter4Reconstructioncoeff(Z_i_mini,theta_inc,phi_inc,nu_max,k_max,k,b,z_explane);
        [Q_S,nu_vect,l_b_vect,Reconcoeff,Z_i_mini]  = targetbeam_interaction(Z_i_mini,b,nu_max,k_max,z_explane,c_1,r_1,Reconcoeff);
        [Gamma_new_i_Q_Mat,A_i_Q,xi_1,xi_2,eta_new_i1,eta_new_i2,sigma_i] =  Gamma_i_rot(Z_i_mini,l_b_vect,b,k_max,nu_max,nu_vect);
        [Gamma_r_mat,eta_new_r1,eta_new_r2,sigma_r] =  GBReflection(xi_1,xi_2,nu_vect,eta_new_i1,eta_new_i2,sigma_i,Gamma_new_i_Q_Mat,l_b_vect);
        [Gamma_new_r_mat,Gamma_new_r,theta_sig_r,phi_sig_r] = N2FF_Gammar_rot(Gamma_r_mat,eta_new_r1,eta_new_r2,sigma_r);
        Z_r = [Q_S,theta_sig_r,phi_sig_r,l_b_vect,A_i_Q,Gamma_new_r]; 
        [sum_sph_mnst_mini, Z_r_mini] = N2FF_GBSM(Gamma_new_r_mat,Gamma_new_r,l_b_vect,theta_obser,phi_obser,phi_sig_r,theta_sig_r,Q_S,k,A_i_Q,Reconcoeff,Z_r);
        Z_r_mini_batches = [Z_r_mini_batches;Z_r_mini];
        Z_i_mini_batches = [Z_i_mini_batches;Z_i_mini];
        sum_sph_mnst = sum_sph_mnst + sum_sph_mnst_mini;
    end
    %% Computation of the relevant number of beams in an observation direction for the lowest frequency in the band 
    eps = 1e-6;
    freq_ind_par  = Z_r_mini_batches(:,12); % Calculates the frequency independent parameters in the system
    Amp_beam = Z_r_mini_batches(:,13);
    [Reconcoeff] = Chapter4Reconstructioncoeff(Z_i_mini_batches,theta_inc,phi_inc,nu_max,k_max,k/2,b,z_explane);
    n_1_2_m_1_2 = Z_i_mini_batches(abs(Reconcoeff.*Amp_beam.*exp(1j*k/2.*freq_ind_par))>eps,1:4);
    n_1 = n_1_2_m_1_2(:,1);
    n_2 = n_1_2_m_1_2(:,2);
    m_1 = n_1_2_m_1_2(:,3);
    m_2 = n_1_2_m_1_2(:,4);
    n_1 = unique(n_1(:).');
    n_2 = unique(n_2(:).');
    m_1 = unique(m_1(:).');
    m_2 = unique(m_2(:).');
    % These parameters will be used to contsruct the beam lattice for the next
    % frequency band. 
    %% Fast RCS in the band calculations:
    k_temp_low = k_max/2; % minimum frequecny of the single octave band
    b_temp = b; % temp value of b, since it changes with the scheme
    k_iter = linspace(k_temp_low,2*k_temp_low,2*del_k); % Single octave
    freq_ind_par  = Z_r_mini_batches(:,12); % Calculates the frequency independent parameters in the system
    Amp_beam = Z_r_mini_batches(:,13);
    [Reconcoeff] = Chapter4Reconstructioncoeff(Z_i_mini_batches,theta_inc,phi_inc,nu_max,k_max,k_iter,b,z_explane);
    RCS_array_oneoct = sum(Reconcoeff.*Amp_beam.*exp(1j*k_iter.*freq_ind_par));
    RCS_date = (Reconcoeff.*Amp_beam.*exp(1j*k_iter.*freq_ind_par));%Delete this afterwards
    [N_arr,beam_array_eps1_mnst_arr,err] = NoBeams_calc(RCS_date);%Delete this afterwards
    %% Computing the next band for the frequency calculations; 
    %upsamp_sch = 1; %X-upsamp  represents 1 and Xi-upsampling represents 
    n_1_min = min(n_1);
    n_1_max = max(n_1);
    n_2_min = min(n_2);
    n_2_max = max(n_2);
    m_1_min = min(m_1);
    m_1_max = max(m_1);
    m_2_min = min(m_2);
    m_2_max = max(m_2);
    if upsamp_sch == 1
        m_1_new = 2*m_1_min:1:2*m_1_max;
        m_2_new = 2*m_2_min:1:2*m_2_max;
        n_1_new = n_1_min:1:n_1_max;
        n_2_new = n_2_min:1:n_2_max;
    else
        m_1_new = m_1_min:1:m_1_max;
        m_2_new = m_2_min:1:m_2_max;
        n_1_new = 2*n_1_min:1:2*n_1_max;
        n_2_new = 2*n_2_min:1:2*n_2_max;
    end
    [N_1,N_2,M_1,M_2] = ndgrid(n_1_new,n_2_new,m_1_new,m_2_new);
    Z_next = [N_1(:),N_2(:),M_1(:),M_2(:)];% Generates all the possible combination of n1,n2,m1,m2
end
