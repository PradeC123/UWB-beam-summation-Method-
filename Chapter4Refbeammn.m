%% Chapter 4 Part 2.5:  %% Reflected beam as a function of m and n 
% This code is dedicated to calculate the RCS of the target at any given
% observation point
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
function Refbeammn = Chapter4Refbeammn(Z,k,k_max,theta_inc,phi_inc,theta,phi,z_b,a,b,nu_max)
    % Finding the contribution at some x = r sin(theta) and z = r cos(theta) 
    %%%%
    x_bar = sqrt(2*pi*nu_max*b/k_max);%Frequency independent lattice
    xi_bar = sqrt(2*pi*nu_max/(b*k_max));
    %%%%
    M1 = floor(2*a/x_bar);
    M2 = floor(2*a/x_bar);
    %M1 = 1;
    %M2 = 1;
    m_1 = -M1:1:M1;
    %m_1 = M1;
    m_2 = -M2:1:M2;%
    %m_2 = M2;
    %n_1max = round(1/xi_bar*sin(theta_inc)*cos(phi_inc));%
    %n_2max = round(1/xi_bar*sin(theta_inc)*sin(phi_inc));
    %C = 64; %Optimization parameter
    %del_n = round(sqrt(C*k_max/(pi*nu*k)));
    %n_1 = n_1max - del_n:1:n_1max + del_n;
    %n_2 = n_2max -del_n:1:n_2max + del_n;
    N1 = floor(sin(60*pi/180)/xi_bar);
    N2 = floor(sin(60*pi/180)/xi_bar);
    %N1 = 1;
    %N2 = 1;
    n_1 = -N1:1:N1;
    n_2 = -N2:1:N2;
    %n_1 = 1;
    [N_1,N_2] = ndgrid(n_1,n_2);
    Z = [N_1(:),N_2(:)];% Generates all the possible combination of n1,n2,m1,m2
    Z(find((Z(:,1)*xi_bar).^2 + (Z(:,2)*xi_bar).^2 >1),:) = [];% Elimination of the beams propagating in evanescent mode
    theta_n = acos(sqrt(1-xi_bar^2*(Z(:,2).^2+Z(:,1).^2)));
    phi_n = mod(atan2(Z(:,2),Z(:,1)),2*pi);
    %%%%
    %%%%
    r = 3*a;
    c1 = [0,0,0];
    x_0 = 0; y_0 = 0; z_0 = 0; 
    x = r*sin(theta)*cos(phi)+c1(1);% Finding the contribution of beac beams at the location.  
    z = r*cos(theta)+c1(3);%
    r1 = a;
    B_contri = zeros(length(m_1),length(n_1));
    A_mu_contri = zeros(length(m_1),length(n_1));
    %c2 = [0,3*a];
    %r2 = a;
    for cnt1 = 1:1:length(m_1)
        for cnt2 = 1:1:length(n_1)
            A_mu_contri(cnt1,cnt2) = nu_max*sqrt(2)*(k/k_max)*exp(-k*b*(sin(theta_inc)*cos(phi_inc)-n_1(cnt2)*xi_bar).^2/2).*...
            exp(-1j*k*sin(theta_inc)*cos(phi_inc)*m_1(cnt1)*x_bar)*exp(-1j*k*z_b*cos(theta_inc)); 
            X_0_mu = [m_1(cnt1)*x_bar,0,z_b];
            theta_i = acos(sqrt(1-xi_bar^2*(n_1(cnt2)^2)));
            phi_i = atan2(0,n_1(cnt2));
            Gamma_i = [1/(1j*b*(cos(theta_i))^2),0;0,1/(1j*b)];
            [QS,theta_r,phi_r,l_b,Gamma_r,A_Q,hit] = multi_intersection_3D(X_0_mu,theta_i,phi_i,Gamma_i,c1,r1);%;    
            if hit == 0 
                B_contri(cnt1,cnt2) = 0; %THe beam doesnot intersect the sphere
            else
                r_xyz = [x - QS(1);   -QS(2); z - QS(3)];
                R_etai1etai2sig = [cos(theta_r)*cos(phi_r),cos(theta_r)*sin(phi_r),-sin(theta_r);...
                -sin(phi_r),cos(phi_r),0; sin(theta_r)*cos(phi_r),sin(theta_r)*sin(phi_r),cos(theta_r)];
                etar1etar2sigr = R_etai1etai2sig* r_xyz;
                eta_r1 = etar1etar2sigr(1);
                eta_r2 = etar1etar2sigr(2);
                sigma_r = etar1etar2sigr(3);
                Gamma_r =  transpose(reshape(Gamma_r,[2,2]));
                eta_r = [eta_r1;eta_r2];
                if (sigma_r >=0.1*r)
                    Gamma_r_sig = inv(inv(Gamma_r) + sigma_r*eye(2));
                    B_contri(cnt1,cnt2) = A_Q*sqrt(det(Gamma_r_sig)/det(Gamma_r))*exp(-1j*k*sigma_r)*exp(-1j*k/2*transpose(eta_r)*Gamma_r_sig*eta_r);
                end
             end
        end
    end 
    Fig1 = figure(1);
    imagesc(n_1*xi_bar,flipud(m_1*x_bar),20*log10(abs(B_contri)));
    set(gca,'YDir','normal');
    axis tight;
    %imagesc(n_1,n_2,a_mu_coeff);
    colormap(jet(256));
    caxis([-140,0]);
    colorbar;
    xlabel('$\xi_n$','Interpreter','latex');
    ylabel("x_m");
    set(gca,'FontSize',14);
    ax = gca; 
    ax.XTickMode = 'manual';
    ax.YTickMode = 'manual';
    ax.ZTickMode = 'manual';
    ax.XLimMode = 'manual';
    ax.YLimMode = 'manual';
    ax.ZLimMode = 'manual';
    saveas(Fig1,'UWB-PS-BSM-Bst-intlcedxixsamp-bmu-j5.png');
    Fig2 = figure(2);
    imagesc(n_1*xi_bar,flipud(m_1*x_bar),20*log10(abs(A_mu_contri)));
    set(gca,'YDir','normal');
    axis tight;
    %imagesc(n_1,n_2,a_mu_coeff);
    colormap(jet(256));
    caxis([-140,0]);
    colorbar;
    xlabel('$\xi_n$','Interpreter','latex');
    ylabel("x_m");
    set(gca,'FontSize',14);
    ax = gca; 
    ax.XTickMode = 'manual';
    ax.YTickMode = 'manual';   
    ax.ZTickMode = 'manual';
    ax.XLimMode = 'manual';
    ax.YLimMode = 'manual';
    ax.ZLimMode = 'manual';
    saveas(Fig2,'UWB-PS-BSM-Bst-intlcedxixsamp-amu-j5.png');
    Fig3 = figure(3);
    imagesc(n_1*xi_bar,flipud(m_1*x_bar),20*log10(abs(A_mu_contri.*B_contri)));
    set(gca,'YDir','normal');
    axis tight;
    %imagesc(n_1,n_2,a_mu_coeff);
    colormap(jet(256));
    caxis([-140,0]);
    colorbar;
    xlabel('$\xi_n$','Interpreter','latex');
    ylabel("x_m");
    set(gca,'FontSize',14);
    ax = gca; 
    ax.XTickMode = 'manual';
    ax.YTickMode = 'manual';
    ax.ZTickMode = 'manual';
    ax.XLimMode = 'manual';
    ax.YLimMode = 'manual';
    ax.ZLimMode = 'manual';
    saveas(Fig3,'UWB-PS-BSM-Bst-intlcedxixsamp-amubmu-j5.png');
end