%% Accuracy as a function of ka and theta 
% This code deals with the error between the absolute value difference
% between the scattered field by GBSM and PO as a functioin of theta for
% several values of ka, usually the method fails at roughly theta = 60 deg
%% Code
function [scatt,RCS] = Errscatt(k,k_max,nu_max,a,b,theta_inc,phi_inc,theta_obser,phi_obser)
    x_bar = sqrt(2*pi*nu_max*b/k_max);%Frequency independent lattice
    xi_bar = sqrt(2*pi*nu_max/(b*k_max));
    z_b = -1*a; %location of the expansion plane 
    %% Lattice expansion indices.
    M1 = floor(1.05*a/x_bar);
    M2 = floor(1.05*a/x_bar);
    m_1 = -M1:1:M1;
    m_2 = -M2:1:M2;%
    N1 = floor(sin(45*pi/180)/xi_bar);
    N2 = floor(sin(45*pi/180)/xi_bar);
    n_1 = -N1:1:N1;
    n_2 = -N2:1:N2;
    [N_1,N_2,M_1,M_2] = ndgrid(n_1,n_2,m_1,m_2);
    Z = [N_1(:),N_2(:),M_1(:),M_2(:)];% Generates all the possible combination of n1,n2,m1,m2
    Z(find((Z(:,1)*xi_bar).^2 + (Z(:,2)*xi_bar).^2 >1),:) = [];% Elimination of the beams propagating in evanescent mode
    theta_n = acos(sqrt(1-xi_bar^2*(Z(:,2).^2+Z(:,1).^2)));
    phi_n = mod(atan2(Z(:,2),Z(:,1)),2*pi);
    Z(:,5) = z_b; 
    Z(:,6) = theta_n;
    Z(:,7) = phi_n; 
    m1 = Z(:,3);
    m2 = Z(:,4);
    n1 = Z(:,1);
    n2 = Z(:,2);
     %% Expansion coefficiens calculation on z = 0 plane
    % Neglecting the edge effects[difraction manifolds]... 
    a_mu = nu_max^2*2*(k/k_max)^2*exp(-k*b*(sin(theta_inc)*cos(phi_inc)-Z(:,1)*xi_bar).^2/2).*exp(-k*b*(sin(theta_inc)*sin(phi_inc)-Z(:,2)*xi_bar).^2/2).*...
        exp(-1j*k*sin(theta_inc)*cos(phi_inc)*Z(:,3)*x_bar).*exp(-1j*k*sin(theta_inc)*sin(phi_inc)*Z(:,4)*x_bar).*exp(-1j*k*z_b*cos(theta_inc)); 
    Z(:,8) = a_mu;
    Z(find(abs(a_mu)<10^-7),:)=[];
    %% Computing the hit point on the sphere.
    tic;
    Z_r = zeros(length(Z(:,1)),11);
    QS = zeros(length(Z(:,1)),3);
    theta_r = zeros(length(Z(:,1)),1);
    phi_r = zeros(length(Z(:,1)),1);
    l_b = zeros(length(Z(:,1)),1);
    hit = zeros(length(Z(:,1)),1);
    Gamma_r = zeros(length(Z(:,1)),4);
    A_Q = zeros(length(Z(:,1)),1);
    c1 = [0,0,0];
    r1 = a;
    for cnt1 = 1:1:length(Z(:,1))
        m_i1 = Z(cnt1,3);
        m_i2 = Z(cnt1,4); 
        X_0 = [m_i1*x_bar,m_i2*x_bar,z_b];
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
    Z_r(:,12) = hit;%(n,m) that hits the sphere
    Z_r(:,13) = Z(:,8);% a_mu>10^-6; Truncation paramter used
    Z_r(find(hit==0),:) = [];% Deleteing all the elements that corresponds to hit = 0
    %% FF Monostatic Calculations 
    Sph_mnst = 0;
    r = 10000000000000*a;
    for cnt1 = 1:1:length(Z_r(:,1)) 
        QS = Z_r(cnt1,1:3);
        theta_r = Z_r(cnt1,4); 
        phi_r = Z_r(cnt1,5);
        l_b = Z_r(cnt1,6); 
        Gamma_r = Z_r(cnt1,7:10);
        Gamma_r =  transpose(reshape(Gamma_r,[2,2]));
        A_Q = Z_r(cnt1,11);
        a_mu = Z_r(cnt1,13);
        r_hat = [sin(theta_obser)*cos(phi_obser),sin(theta_obser)*sin(phi_obser),cos(theta_obser)]; % Monostatic observation angle
        theta = acos(r_hat(3));
        phi = mod(atan2(r_hat(2),r_hat(1)),2*pi);
        sig_dir = sin(theta_r)*cos(phi_r)*(r*sin(theta)*cos(phi)-QS(1)) + sin(theta_r)*sin(phi_r)*(r*sin(theta)*sin(phi)-QS(2))...
            + cos(theta_r)*(r*cos(theta)-QS(3));
        if (sig_dir>=0.90*r) %In the lit region
            omega = [sin(theta_obser)*cos(theta_r)*cos(phi_obser-phi_r)-cos(theta_obser)*sin(theta_r);sin(theta_obser)*sin(phi_obser-phi_r)];
            ff_GB_decay = exp(1j*k/2*transpose(omega)*((Gamma_r)^-1)*omega);% Far field Gaussian Beam angular decay term
            ff_phase_shift = exp(1j*k*dot(r_hat,QS)); % Far field fraunhofer phase shift term
            Psi_thet_phi_mnst = A_Q*sqrt(1/(det(Gamma_r)))*ff_GB_decay*ff_phase_shift;
            Sph_mnst = Sph_mnst + a_mu*Psi_thet_phi_mnst*exp(-1j*k*l_b);
            Z_r(cnt1,14) = a_mu*Psi_thet_phi_mnst*exp(-1j*k*l_b);
        end
    end
    %% Results
    scatt = Sph_mnst;
    RCS = 4/(a^2)*abs(Sph_mnst).^2;   
end