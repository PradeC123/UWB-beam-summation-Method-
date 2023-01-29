%% Expansion parameters
tic
a = 1;
b = 1;
k_max = 8000;
nu_max = 0.33;
z_explane = 0;
k = 8000;
theta_inc = 0;
phi_inc = 0;
%% Parameters of the sphere
c_1 = [0,0,0]; %The location of the center of the sphere 
r_1 = 1; % The radius of the sphere 
Tar_surf = [c_1,r_1];%Target surface
theta_obser = pi; %Observation angles in the far zone
phi_obser = pi;
%% RCS calculation for one single frequency band 
Expplane = Chapter4Expansionplane(a,b,k_max,nu_max,z_explane,theta_inc,phi_inc);
Z = Expplane;
del = 2*1e7;
mini_batch = floor(length(Expplane)/del);
sum_sph_mnst = 0; 
Z_i_mini_batches = [];
Z_r_mini_batches = [];
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
    Z_i_mini = Z(i_cnt*del+1:length(Expplane),:);
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
    disp(4*abs(sum_sph_mnst).^2)
end
toc 
%% Computation of the relevant number of beams in an observation direction for the lowest frequency in the band 
tic;
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
toc;
%% Fast RCS in the band calculations: 
tic;
k_min = 1000;
k_temp_low = k_min; % minimum frequecny of the band
k_temp_max = k_max; % minimum frequecny of the band
b_temp = b; % temp value of b, since it changes with the scheme
i = 1; 
del_k = 125/2^3;
k_iter = k_temp_low:2^i*del_k:k_temp_max; % Single octave
freq_ind_par  = Z_r_mini_batches(:,12); % Calculates the frequency independent parameters in the system
Amp_beam = Z_r_mini_batches(:,13);
[Reconcoeff] = Chapter4Reconstructioncoeff(Z_i_mini_batches,theta_inc,phi_inc,nu_max,k_max,k_iter,b,z_explane);
sum(Reconcoeff.*Amp_beam.*exp(1j*k_iter.*freq_ind_par));
toc;
%% Computing the next band for the frequency calculations; 
x_upsamp = 1;
n_1_min = min(n_1);
n_1_max = max(n_1);
n_2_min = min(n_2);
n_2_max = max(n_2);
m_1_min = min(m_1);
m_1_max = max(m_1);
m_2_min = min(m_2);
m_2_max = max(m_2);
if x_upsamp == 1
    m_1_new = 2*m_1_min:1:2*m_1_max;
    m_2_new = 2*m_2_min:1:2*m_2_max;
    n_1_new = n_1_min:1:n_1_max;
    n_2_new = n_2_min:1:n_2_max;
    length(m_1_new)^2*length(n_1_new)^2;
else
    m_1_new = m_1_min:1:m_1_max;
    m_2_new = m_2_min:1:m_2_max;
    n_1_new = 2*n_1_min:1:2*n_1_max;
    n_2_new = 2*n_2_min:1:2*n_2_max;
    length(m_1_new)^2*length(n_1_new)^2;
end
[N_1,N_2,M_1,M_2] = ndgrid(n_1,n_2,m_1,m_2);
Z_next = [N_1(:),N_2(:),M_1(:),M_2(:)];% Generates all the possible combination of n1,n2,m1,m2
Z_next(find((Z_next(:,1)*xi_bar).^2 + (Z_next(:,2)*xi_bar).^2 >1),:) = [];% Elimination of the beams propagating in evanescent mode
theta_n = acos(sqrt(1-xi_bar^2*(Z_next(:,2).^2+Z_next(:,1).^2)));
phi_n = mod(atan2(Z_next(:,2),Z_next(:,1)),2*pi);
Z_next(:,5) = z_b; 
Z_next(:,6) = theta_n;
Z_next(:,7) = phi_n; 
%m1 = Z_next(:,3);
%m2 = Z_next(:,4);
%n1 = Z_next(:,1);
%n2 = Z_next(:,2);
%Expplane = Z_next;
%% Multiband analysis. 
[Z_next,RCS_array_oneoct] = beam_RCS_oneoct(Z,a,b,k_max,nu_max,z_explane,theta_inc,phi_inc,theta_obser,phi_obser,del_k,upsamp_sch); 
k_min = 1000;
k_max = 16000;
del_k = 125/8;
a = 1; 
b_bar = a;
z_b = 0;
J = log2(k_max/k_min); % Total number of subbands in the entire frequency band.
RCS_array_multiband = [];
k_subbands = []; 
for i = 1:1:J
    k_max = 2^i*k_min;
    k_temp_low = k_max/2;
    k_iter = k_temp_low:2^i*del_k:2*k_temp_low;
    if i == 1
        tic
        Z = Chapter4Expansionplane(a,b_bar,k_max,nu_max,z_explane,theta_inc,phi_inc);
        [Z_next,RCS_array_oneoct] = beam_RCS_oneoct(Z,a,b_bar,k_max,nu_max,z_explane,theta_inc,phi_inc,theta_obser,phi_obser,del_k,upsamp_sch);
        k_subbands = [k_subbands,k_iter];
        RCS_array_multiband = [RCS_array_multiband,RCS_array_oneoct];
        toc 
    end
    if i > 1 %Beam expansion for the next lattice
        xi_bar = sqrt(2*pi*nu_max/(b_temp*k_max));
        Z_next(find((Z_next(:,1)*xi_bar).^2 + (Z_next(:,2)*xi_bar).^2 >1),:) = [];% Elimination of the beams propagating in evanescent mode
        theta_n = acos(sqrt(1-xi_bar^2*(Z_next(:,2).^2+Z_next(:,1).^2)));
        phi_n = mod(atan2(Z_next(:,2),Z_next(:,1)),2*pi);
        Z_next(:,5) = z_b;
        Z_next(:,6) = theta_n;
        Z_next(:,7) = phi_n;
        Z_next(find((Z_next(:,1)*xi_bar).^2 + (Z_next(:,2)*xi_bar).^2 >1),:) = [];% Elimination of the beams propagating in evanescent mode
        theta_n = acos(sqrt(1-xi_bar^2*(Z_next(:,2).^2+Z_next(:,1).^2)));
        phi_n = mod(atan2(Z_next(:,2),Z_next(:,1)),2*pi);
        Z_next(:,5) = z_b;
        Z_next(:,6) = theta_n;
        Z_next(:,7) = phi_n;
        [Z_next,RCS_array_oneoct] = beam_RCS_oneoct(Z,a,b_temp,k_max,nu_max,z_explane,theta_inc,phi_inc,theta_obser,phi_obser,del_k,upsamp_sch);
        k_subbands = [k_subbands,k_iter];
        RCS_array_multiband = [RCS_array_multiband,RCS_array_oneoct];
    end
    %Beam expansion for the next lattice 
    switch swtch_scheme %Switch cas es  to compute the RCS of the target by all 4 schemes successfully deployed 
     case 1
         b_temp = b_bar*(2^(-i));% X-sampling scheme
         upsamp_sch = 1; 
     case 2
         b_temp = b_bar*(2^(i));% X-sampling scheme
         upsamp_sch = 0; 
     case 3
         if mod(i+1,2) == 0
             b_temp = b_bar/2;
             upsamp_sch = 1;
         else
             b_temp = b_bar;
             upsamp_sch = 0;
         end
     case 4
         if mod(i+1,2) == 0
             b_temp = 2*b_bar;
             upsamp_sch = 0; 
         else
             b_temp = b_bar;
             upsamp_sch = 1; % x upsampling
         end
    end
end
%%
fig1 = figure(1);
yyaxis left;
plot(k_iter, 4/a^2*abs(sum(Reconcoeff.*Amp_beam.*exp(1j*k_iter.*freq_ind_par))).^2,'LineWidth',2);
ylim([0.9 1.1]);
yyaxis right;
%plot(k_iter, 20*log10(abs( exp(1j*2*k_iter)/2 - sum(Reconcoeff.*Amp_beam.*exp(1j*k_iter.*freq_ind_par)))),'-.','LineWidth',2);
plot(k_iter, 20*log10(abs(1-4/a^2*abs(sum(Reconcoeff.*Amp_beam.*exp(1j*k_iter.*freq_ind_par))).^2)),'-.','LineWidth',2);
ylim([-80 -10]);
grid on;
set(gca,'FontSize',14);
yyaxis left;
ylabel('$\frac{\sigma}{\pi a^2}$','Interpreter','latex','FontSize',20);
xlim([1000,8000])
xlabel('{\itk}');
saveas(fig1,'Chapter4_Sphsim_nuk1000.png');



