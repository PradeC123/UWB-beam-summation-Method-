%% Initial parameters 
%clear all;
%clc;
theta_inc = pi/9; %inicident angle of propagation of the field 
phi_inc = 0; %incident angle of propagation of the field 
i = 4;
upsamp = 4;
k_max = 2000*2^(i-1); % The normalized wavenumber also the highest frequency in the frequency band
%k = 1500*2^(i-1);  %k = 0.5 0.75 1 observation frequency
k = sqrt(1000*2^(i-1)*2000*2^(i-1));
a = 1;% The radius of the sphere 
b = (upsamp==1)*a*2^(-i+1) + (upsamp==2)*a*2^(i-1) + (upsamp==3)*a*2^(-1*mod(i-1,2)) + (upsamp==4)*a*2^(mod(i-1,2));
W = sqrt(b/k);
nu_max = 0.20;
x_bar = sqrt(2*pi*nu_max*b/k_max);%Frequency independent lattice
xi_bar = sqrt(2*pi*nu_max/(b*k_max));
Theta_d1 = sqrt(1/(k*b));
z_b = -1*a; %location of the expansion plane 
%% Lattice expansion indices.
%L = x1_max; %the width of the square aperture in x1 axis. Onl relevant for truncation of the beams 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Respective indices in m_1,m_2,n_1,n_2
M1 = floor(1*a/x_bar);
M2 = floor(1*a/x_bar);
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
%n_2 = 1;
%% Reflected beam as a function of m1 and n1
% Finding the contribution at some x = r sin(theta) and z = r cos(theta)
phi_inc = 0; 
theta_inc = pi/9;
%theta = pi/180*180-theta_inc;
theta = 5*pi/6;
phi = pi;
r = 3*a;
c1 = [0,0,0];
x_0 = 0; y_0 = 0; z_0 = 0; 
x = r*sin(theta)*cos(phi)+c1(1);% Finding the contribution of beam beams at the location.  
y = r*sin(theta)*sin(phi)+c1(2);
z = r*cos(theta)+c1(3);%
r1 = a;
B_contri = zeros(length(m_1),length(n_1));
A_mu_contri = zeros(length(m_1),length(n_1));
%c2 = [0,3*a];
%r2 = a;
m_2_max = 0; 
n_2_max = 0; 
for cnt1 = 1:1:length(m_1)
    for cnt2 = 1:1:length(n_1)
        A_mu_contri(cnt1,cnt2) = nu_max*sqrt(2)*(k/k_max)*exp(-k*b*(sin(theta_inc)*cos(phi_inc)-n_1(cnt2)*xi_bar).^2/2).*...
        exp(-1j*k*sin(theta_inc)*cos(phi_inc)*m_1(cnt1)*x_bar)*exp(-1j*k*z_b*cos(theta_inc)); 
        X_0_mu = [m_1(cnt1)*x_bar,m_2_max*x_bar,z_b];
        theta_i = acos(sqrt(1-xi_bar^2*(n_1(cnt2)^2+n_2_max^2)));
        phi_i = atan2(n_2_max,n_1(cnt2));
        Gamma_i = [1/(1j*b*(cos(theta_i))^2),0;0,1/(1j*b)];
        [QS,theta_r,phi_r,l_b,Gamma_r,A_Q,hit] = multi_intersection_3D(X_0_mu,theta_i,phi_i,Gamma_i,c1,r1);%;    
        if hit == 0 
            B_contri(cnt1,cnt2) = 0; %The beam doesnot intersect the sphere
        else
            r_xyz = [x - QS(1); y - QS(2); z - QS(3)];
            R_etai1etai2sig = [cos(theta_r)*cos(phi_r),cos(theta_r)*sin(phi_r),-sin(theta_r);...
            -sin(phi_r),cos(phi_r),0; sin(theta_r)*cos(phi_r),sin(theta_r)*sin(phi_r),cos(theta_r)];
            etar1etar2sigr = R_etai1etai2sig* r_xyz;
            eta_r1 = etar1etar2sigr(1);
            eta_r2 = etar1etar2sigr(2);
            sigma_r = etar1etar2sigr(3);
            Gamma_r =  transpose(reshape(Gamma_r,[2,2]));
            eta_r = [eta_r1;eta_r2];
            if (sigma_r >=0)
                Gamma_r_sig = inv(inv(Gamma_r) + sigma_r*eye(2));
                B_contri(cnt1,cnt2) = A_Q*sqrt(det(Gamma_r_sig)/det(Gamma_r))*exp(-1j*k*sigma_r)*exp(-1j*k/2*transpose(eta_r)*Gamma_r_sig*eta_r);
            end
         end
    end
end 
%%
Fig3 = figure(3);
imagesc(n_1*xi_bar,flipud(m_1)*x_bar,20*log10(abs(A_mu_contri.*B_contri)));
set(gca,'YDir','normal');
axis tight;
%imagesc(n_1,n_2,a_mu_coeff);
colormap(jet(256));
caxis([-120,0]);
colorbar;
xlabel('n1');
ylabel('m1');
ax = gca; 
ylabel('$x_m$','Interpreter','latex','FontSize',14);
xlabel('$\xi_n$','Interpreter','latex','FontSize',14);
%% 
B_contr = (abs(A_mu_contri.*B_contri))>=10^-6;
%B_contr = (abs(A_mu_contri.*B_contri));
M_tot = sum((sum(B_contr,2)>0)==1);
N_tot = sum((sum(B_contr,1)>0)==1);
n_1_arr = n_1(find((sum(B_contr,1)>0)==1))*xi_bar;
m_1_arr = m_1(find((sum(B_contr,2)>0)==1))*x_bar;
hold on;
rectangle('Position',[n_1_arr(1),m_1_arr(1),N_tot*xi_bar,M_tot*x_bar],'EdgeColor','b','LineWidth',3);
xlim([0.2,0.5]);
ylim([-0.8,-0.2]);
M_tot*x_bar;
N_tot*xi_bar;
