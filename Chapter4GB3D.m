% Date - 02/09/2020 (18:42) Israel Daylight Time (Started on)
%% Initial parameters 
%clear all;
%clc;
theta_inc = 0; %inicident angle of propagation of the field 
phi_inc = 0; %incident angle of propagation of the field 
k_max = 000; % The normalized wavenumber also the highest frequency in the frequency band
k = 12000;  %k = 0.5 0.75 1 observation frequency
a = 1;% The radius of the sphere 
b = 8*a; %b = 1000,2000,4000, Collimation distance in the frequency normalized system
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
%n_2 = 1;
%%
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
% a_mu = Rec_coeff_PW(theta_inc,phi_inc,nu,b,k,m1,m2,n1,n2,z_b);
Z(:,8) = a_mu;
%Z(find((abs(a_mu))<10^-15),:)=[]; 
Z(find(abs(a_mu)<10^-10),:)=[];
%% PLotting the expansion coefficient
 a_mu_coeff = nu_max^2*2*(k/k_max)^2*exp(-k*b*(sin(theta_inc)*sin(phi_inc)-transpose(n_2)*xi_bar).^2/2)...
    .*exp(-k*b*(sin(theta_inc)*cos(phi_inc)-(n_1)*xi_bar).^2/2); %plotting the expansion coefficents as a fucntion of n_1*n_2 
Fig0 = figure(1);
plot(0);
hold on
imagesc(n_1,n_2,20*log10(a_mu_coeff/max(max(a_mu_coeff))));
axis tight;
%imagesc(n_1,n_2,a_mu_coeff);
colormap(jet(256));
caxis([-140,0]);
colorbar;
xlabel("n_1");
ylabel("n_2");
set(gca,'FontSize',10);
%% Incident field as a function of m1 and n1
% Finding the contribution at some x = r sin(theta) and z = r cos(theta)
phi_inc = 0; 
theta_inc = pi/12;
x = x_bar/2;% Finding the contribution of beam beams at the location.  
y = 0;
z = b/2;%
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
        r_xyz = [x- m_1(cnt1)*x_bar; y-m_2_max*x_bar; z - z_b];
        R_etai1etai2sig = [cos(theta_i)*cos(phi_i),cos(theta_i)*sin(phi_i),-sin(theta_i);...
            -sin(phi_i),cos(phi_i),0; sin(theta_i)*cos(phi_i),sin(theta_i)*sin(phi_i),cos(theta_i)];
        r_etasig = R_etai1etai2sig*r_xyz;
        eta  = [r_etasig(1),r_etasig(2)];
        sig_i = r_etasig(3);
        Gamma_sig_i = inv(inv(Gamma_i)+sig_i*eye(2));
        Amp_EC = sqrt(det(Gamma_sig_i)/det(Gamma_i)); %conservation of energy
        lin_phase = exp(-1j*k*sig_i); %linear phase term
        para_phase_amp = exp(-1j*k/2*eta(1)^2*Gamma_sig_i(1,1))*exp(-1j*k/2*eta(2)^2*Gamma_sig_i(2,2)); %paraxial phase and amplitude correction term for the Gaussian beam
        B_contri(cnt1,cnt2) = Amp_EC*lin_phase*para_phase_amp;
    end
end
Fig1 = figure(1);
imagesc(n_1*xi_bar,flipud(m_1*x_bar),20*log10(abs(B_contri)));
set(gca,'YDir','normal');
axis tight;
colormap(jet(256));
caxis([-120,0]);
colorbar;
xlabel('n1');
ylabel('m1');
ylabel('$x_m$','Interpreter','latex','FontSize',14);
xlabel('$\xi_n$','Interpreter','latex','FontSize',14);% ylim([-15,15]);
%saveas(Fig1,'monostatic-bmu-j1-2m1n1.png');
ax = gca; 
Fig2 = figure(2);
imagesc(n_1*xi_bar,flipud(m_1*x_bar),20*log10(abs(A_mu_contri)));
set(gca,'YDir','normal');
axis tight;
colormap(jet(256));
caxis([-120,0]);
colorbar;
ylabel('$x_m$','Interpreter','latex','FontSize',14);
xlabel('$\xi_n$','Interpreter','latex','FontSize',14);
%saveas(Fig2,'monostatic-amu-j1-2m1n1.png');
Fig3 = figure(3);
imagesc(n_1,flipud(m_1*x_bar),20*log10(abs(A_mu_contri.*B_contri)));
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
%saveas(Fig3,'monostatic-ambmu-j1-2m1n1.png');
%%
% In this section we are ploting the phase for x = -a:dx:a, y=0 and z=0.
% long along the x axis along the aparture plane.
x = -x_bar:1:x_bar;% discretization of 0.1*lambda.
%x = 0;
y = 0; z= a;% The case to analyse x = -a to x = a and z = 0, initially.
U_tot = zeros(1,length(x));
U_tot1 = zeros(1,length(x));
sum1 = 0;
U_tot2 = zeros(1,length(x));
sum2 = 0;
U_tot3 = zeros(1,length(x));
U_tot4 = zeros(1,length(x));
U_tot5 = zeros(1,length(x));
U_tot6 = zeros(1,length(x));
sum3 = 0;
lcl_field = zeros(1,length(x));
for cnt2 = 1:1:length(x)
    sum3 = 0;
    for cnt1 = 1:1:length(Z(:,1))
        theta_ni = Z(cnt1,6);
        phi_ni = Z(cnt1,7); 
        m_i1 = Z(cnt1,3); 
        m_i2 = Z(cnt1,4); 
        r_xyz = [x(cnt2) - m_i1*x_bar; y-m_i2*x_bar; z - z_b];
        R_etai1etai2sig = [cos(theta_ni)*cos(phi_ni),cos(theta_ni)*sin(phi_ni),-sin(theta_ni);...
            -sin(phi_ni),cos(phi_ni),0; sin(theta_ni)*cos(phi_ni),sin(theta_ni)*sin(phi_ni),cos(theta_ni)];
        r_etasig = R_etai1etai2sig*r_xyz;
        eta  = [r_etasig(1),r_etasig(2)];
        sig_i = r_etasig(3);
        Gamma_i = [1/(1j*b*cos(theta_ni)^2),0;0,1/(1j*b)];
        Gamma_sig_i = inv(inv(Gamma_i)+sig_i*eye(2));
        Amp_EC = sqrt(det(Gamma_sig_i)/det(Gamma_i)); %conservation of energy
        lin_phase = exp(-1j*k*sig_i); %linear phase term
        para_phase_amp = exp(-1j*k/2*eta(1)^2*Gamma_sig_i(1,1))*exp(-1j*k/2*eta(2)^2*Gamma_sig_i(2,2)); %paraxial phase and amplitude correction term for the Gaussian beam
        B_mu = Amp_EC*lin_phase*para_phase_amp;
        Z(cnt1,9) = B_mu;
        sum3 = sum3 + Z(cnt1,8)*Z(cnt1,9);
    end
    U_tot(cnt2) = sum(Z(:,8).*Z(:,9));
    eps1 = find(abs(Z(:,8).*Z(:,9))>0.1);
    U_tot1(cnt2) = sum(Z(eps1,8).*Z(eps1,9));
    eps2 = find(abs(Z(:,8).*Z(:,9))>0.01);
    U_tot2(cnt2) = sum(Z(eps2,8).*Z(eps2,9));
    eps3 = find(abs(Z(:,8).*Z(:,9))>0.001);
    U_tot3(cnt2) = sum(Z(eps3,8).*Z(eps3,9));
    eps4 = find(abs(Z(:,8).*Z(:,9))>0.0001);
    U_tot4(cnt2) = sum(Z(eps4,8).*Z(eps4,9));
    eps5 = find(abs(Z(:,8).*Z(:,9))>0.00001);
    U_tot5(cnt2) = sum(Z(eps5,8).*Z(eps5,9));
    eps6 = find(abs(Z(:,8).*Z(:,9))>0.000001);
    U_tot6(cnt2) = sum(Z(eps6,8).*Z(eps6,9));
    disp(cnt2);
end
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
% %c2 = [-1*a,0,0];
% %r2 = a; 
% [X1,Y1,Z1] = sphere;
% % r = 200;
% % X2 = X1 * r;
% % Y2 = Y1 * r;
% % Z2 = Z1 * r;
% % surf(X2,Y2,Z2);
% hold on;
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
Z_r(:,12) = hit;%(n,m) that hi ts the sphere
Z_r(:,13) = Z(:,8);% a_mu>10^-6; Truncation paramter used
Z_r(find(hit==0),:) = [];% Deleteing all the elements that corresponds to hit = 0
%Z_r(find(Z_r(:,13)<10^-10),:)=[];
%Z_r(find((abs(Z_r(:,13))/max(abs(Z_r(:,13)))10^-10)),:)=[];
toc;
%% Reflected beam as a function of m1 and n1
% Finding the contribution at some x = r sin(theta) and z = r cos(theta)
phi_inc = 0; 
theta_inc = pi/6;
theta = pi/180*180-theta_inc;
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
            if (sigma_r >=0.5*r)
                Gamma_r_sig = inv(inv(Gamma_r) + sigma_r*eye(2));
                B_contri(cnt1,cnt2) = A_Q*sqrt(det(Gamma_r_sig)/det(Gamma_r))*exp(-1j*k*sigma_r)*exp(-1j*k/2*transpose(eta_r)*Gamma_r_sig*eta_r);
            end
         end
    end
end 
Fig1 = figure(1);
imagesc(n_1*xi_bar,flipud(m_1)*x_bar,20*log10(abs(B_contri)));
set(gca,'YDir','normal');
axis tight;
%imagesc(n_1,n_2,a_mu_coeff);
colormap(jet(256));
caxis([-120,0]);
colorbar;
xlabel('n1');
ylabel('m1');
%set(gca,'visible','off')
colorbar('off')
% ylim([-15,15]);
% xlim([-15,15]);
xlabel('$\xi_n$','Interpreter','latex');
ylabel("x_m");
set(gca,'FontSize',14);
ax = gca; 
% ax.XTickMode = 'manual';
% ax.YTickMode = 'manual';
% ax.ZTickMode = 'manual';
% ax.XLimMode = 'manual';
% ax.YLimMode = 'manual';
% ax.ZLimMode = 'manual';
%saveas(Fig1,'monostatic-bmu-j1-m1n1.png');
Fig2 = figure(2);
imagesc(n_1*xi_bar,flipud(m_1)*x_bar,20*log10(abs(A_mu_contri)));
set(gca,'YDir','normal');
axis tight;
%imagesc(n_1,n_2,a_mu_coeff);
colormap(jet(256));
caxis([-120,0]);
colorbar;
%xlabel('n1');
%ylabel('m1');
%set(gca,'visible','off')
%colorbar('off')
% ylim([-15,15]);
% xlim([-15,15]);
xlabel('$\xi_n$','Interpreter','latex');
ylabel("x_m");
set(gca,'FontSize',14);
% ax = gca; 
% ax.XTickMode = 'manual';
% ax.YTickMode = 'manual';
% ax.ZTickMode = 'manual';
% ax.XLimMode = 'manual';
% ax.YLimMode = 'manual';
% ax.ZLimMode = 'manual';
%saveas(Fig2,'monostatic-amu-j2-m1n1.png');
Fig3 = figure(3);
imagesc(n_1*xi_bar,flipud(m_1)*x_bar,20*log10(abs(A_mu_contri.*B_contri)));
set(gca,'YDir','normal');
axis tight;
%imagesc(n_1,n_2,a_mu_coeff);
colormap(jet(256));
caxis([-120,0]);
colorbar;
%xlabel('n1');
%ylabel('m1');
%colorbar('off')
% ylim([-15,15]);
% xlim([-15,15]);
xlabel('$\xi_n$','Interpreter','latex');
ylabel("x_m");
set(gca,'FontSize',14);
ax = gca; 
%set(gca,'visible','off')
%colorbar('off')
% ax.XTickMode = 'manual';
% ax.YTickMode = 'manual';
% ax.ZTickMode = 'manual';
% ax.XLimMode = 'manual';
% ax.YLimMode = 'manual';
% ax.ZLimMode = 'manual';
%xlim([5,20]);
%ylim([-30 10]);
%saveas(Fig3,'monostatic-ambmu-j2-m1n1.png');
%% Reflected beam as a function of m2 and n2
% Finding the contribution at some x = r sin(theta) and z = r cos(theta)
phi_inc = 0; 
theta_inc = pi/9;
theta = pi/180*180-theta_inc;
phi = pi+phi_inc;
r = 3*a;
c1 = [0,0,0];
x_0 = 0; y_0 = 0; z_0 = 0; 
x = r*sin(theta)*cos(phi)+c1(1);% Finding the contribution of beam beams at the location.  
y = r*sin(theta)*sin(phi)+c1(2);
z = r*cos(theta)+c1(3);%
r1 = a;
B_contri = zeros(length(m_2),length(n_2));
A_mu_contri = zeros(length(m_2),length(n_2));  
%c2 = [0,3*a];
%r2 = a; 
m_1_max = -16;
n_1_max = 16;
for cnt1 = 1:1:length(m_2)
    for cnt2 = 1:1:length(n_2)
        A_mu_contri(cnt1,cnt2) = nu_max*sqrt(2)*(k/k_max)*exp(-k*b*(sin(theta_inc)*sin(phi_inc)-n_2(cnt2)*xi_bar).^2/2).*...
        exp(-1j*k*sin(theta_inc)*sin(phi_inc)*m_2(cnt1)*x_bar)*exp(-1j*k*z_b*cos(theta_inc)); 
        X_0_mu = [m_1_max*x_bar,m_2(cnt1)*x_bar,z_b];
        theta_i = acos(sqrt(1-xi_bar^2*(n_2(cnt2)^2+n_1_max^2)));
        phi_i = atan2(n_2(cnt2),n_1_max);
        Gamma_i = [1/(1j*b*(cos(theta_i))^2),0;0,1/(1j*b)];
        [QS,theta_r,phi_r,l_b,Gamma_r,A_Q,hit] = multi_intersection_3D(X_0_mu,theta_i,phi_i,Gamma_i,c1,r1);%;    
        if hit == 0 
            B_contri(cnt1,cnt2) = 0; %THe beam doesnot intersect the sphere
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
            if (sigma_r >=0.1*r)
                Gamma_r_sig = inv(inv(Gamma_r) + sigma_r*eye(2));
                B_contri(cnt1,cnt2) = A_Q*sqrt(det(Gamma_r_sig)/det(Gamma_r))*exp(-1j*k*sigma_r)*exp(-1j*k/2*transpose(eta_r)*Gamma_r_sig*eta_r);
            end
         end
    end
end 
Fig1 = figure(1);
imagesc(n_2*xi_bar,flipud(m_2)*x_bar,20*log10(abs(B_contri)));
set(gca,'YDir','normal');
axis tight;
%imagesc(n_1,n_2,a_mu_coeff);
colormap(jet(256));
caxis([-120,0]);
% set(gca,'visible','off')
colorbar('off')
% colorbar;
% xlabel('n2');
% ylabel('m2');
% ylim([-15,15]);
% xlim([-15,15]);
%xlabel('$\xi_n$','Interpreter','latex');
%ylabel("x_m");
%set(gca,'FontSize',14);
ax = gca; 
% ax.XTickMode = 'manual';
% ax.YTickMode = 'manual';
% ax.ZTickMode = 'manual';
% ax.XLimMode = 'manual';
% ax.YLimMode = 'manual';
% ax.ZLimMode = 'manual';
% saveas(Fig1,'monostatic-bmu-j1-m2n2.png');
Fig2 = figure(2);
imagesc(n_2*xi_bar,flipud(m_2)*x_bar,20*log10(abs(A_mu_contri)));
set(gca,'YDir','normal');
axis tight;
%imagesc(n_1,n_2,a_mu_coeff);
colormap(jet(256));
caxis([-120,0]);
% set(gca,'visible','off')
colorbar('off')
% colorbar;
% xlabel('n2');
% ylabel('m2');
% ylim([-15,15]);
% xlim([-15,15]);
%xlabel('$\xi_n$','Interpreter','latex');
%ylabel("x_m");
%set(gca,'FontSize',14);
% ax = gca; 
% ax.XTickMode = 'manual';
% ax.YTickMode = 'manual';
% ax.ZTickMode = 'manual';
% ax.XLimMode = 'manual';
% ax.YLimMode = 'manual';
% ax.ZLimMode = 'manual';
% saveas(Fig2,'monostatic-amu-j2-m2n2.png');
Fig3 = figure(3);
imagesc(n_2*xi_bar,flipud(m_2)*x_bar,20*log10(abs(A_mu_contri.*B_contri)));
set(gca,'YDir','normal');
axis tight;
%imagesc(n_1,n_2,a_mu_coeff);
% colormap(jet(256));
caxis([-120,0]);
colormap(jet(256));
set(gca,'visible','off');
colorbar('off')
% colorbar;
% xlabel('n2');
% ylabel('m2');
% ylim([-15,15]);
% xlim([-15,15]);
%xlabel('$\xi_n$','Interpreter','latex');
%ylabel("x_m");
%set(gca,'FontSize',14);
ax = gca; 
% ax.XTickMode = 'manual';
% ax.YTickMode = 'manual';
% ax.ZTickMode = 'manual';
% ax.XLimMode = 'manual';
% ax.YLimMode = 'manual';
% ax.ZLimMode = 'manual';
%saveas(Fig3,'monostatic-ambmu-j2-m2n2.png');
 %% FF Monostatic Calculations 
tic;
Sph_mnst = 0;
Sph_mnst_1 = 0;
theta_obser = 180*pi/180 - theta_inc;
phi_obser = pi + phi_inc;
r = 10000000000000*a;
for cnt1 = 1:1:length(Z_r(:,1)) 
    QS = Z_r(cnt1,1:3);
    theta_r = Z_r(cnt1,4); 
    phi_r = Z_r(cnt1,5);
    %phi_r = phi_r*(phi_r<pi)+(phi_r-2*pi)*(phi_r>pi);
    l_b = Z_r(cnt1,6); 
    Gamma_r = Z_r(cnt1,7:10);
    Gamma_r =  transpose(reshape(Gamma_r,[2,2]));
    A_Q = Z_r(cnt1,11);
    a_mu = Z_r(cnt1,13);
    r_hat = [sin(theta_obser)*cos(phi_obser),sin(theta_obser)*sin(phi_obser),cos(theta_obser)]; % Monostatic observation angle
    theta = acos(r_hat(3));
    phi = mod(atan2(r_hat(2),r_hat(1)),2*pi);
   % phi = atan2(r_hat(2),r_hat(1));
    %phi = phi*(phi>=0) + (2*pi+phi)*(phi<0);
    sig_dir = sin(theta_r)*cos(phi_r)*(r*sin(theta)*cos(phi)-QS(1)) + sin(theta_r)*sin(phi_r)*(r*sin(theta)*sin(phi)-QS(2))...
        + cos(theta_r)*(r*cos(theta)-QS(3));
    if (sig_dir>=0.9*r) %In the lit region
        omega = [sin(theta_obser)*cos(theta_r)*cos(phi_obser-phi_r)-cos(theta_obser)*sin(theta_r);sin(theta_obser)*sin(phi_obser-phi_r)];
        %omega_1 = [omega_1 ; transpose(omega)];
        %omega = [(theta-theta_r); sin(theta_r)*sin(phi-phi_r)];
        %omega_2 = [omega_2 ; transpose(omega)];
        ff_GB_decay = exp(1j*k/2*transpose(omega)*((Gamma_r)^-1)*omega);% Far field Gaussian Beam angular decay term
        ff_phase_shift = exp(1j*k*dot(r_hat,QS)); % Far field fraunhofer phase shift term
        Psi_thet_phi_mnst = A_Q*sqrt(1/(det(Gamma_r)))*ff_GB_decay*ff_phase_shift;
        Sph_mnst = Sph_mnst + a_mu*Psi_thet_phi_mnst*exp(-1j*k*l_b);
        Z_r(cnt1,14) = a_mu*Psi_thet_phi_mnst*exp(-1j*k*l_b);
    end
end
toc
%% Calculating the Bistatic RCS in vector parallelisation manner(Efficient and Fast Method of Computation)
 tic;
r = 100000*a;
Sph_mnst = 0;
theta_obser_arr = transpose(0:0.0005:pi);
phi_obser = pi;
u_bistat_obser_lit = zeros(1,length(theta_obser_arr));
del = 0.80;
for cnt1 = 1:1:length(Z_r(:,1))
    QS = Z_r(cnt1,1:3);
    QS_arr = [transpose(repelem(QS(1),length(theta_obser_arr))),transpose(repelem(QS(2),length(theta_obser_arr)))...
       ,transpose(repelem(QS(3),length(theta_obser_arr)))];
    %QS_arr = [transpose(repelem(QS(1),length(phi_obser))),transpose(repelem(QS(2),length(phi_obser)))...
       % ,transpose(repelem(QS(3),length(phi_obser)))];
    theta_r = Z_r(cnt1,4); 
    phi_r = Z_r(cnt1,5);
    l_b = Z_r(cnt1,6); 
    Gamma_r = Z_r(cnt1,7:10);
    A_Q = Z_r(cnt1,11);
    a_mu = Z_r(cnt1,13);
    r_hat = [sin(theta_obser_arr).*cos(phi_obser),sin(theta_obser_arr).*sin(phi_obser),repelem(cos(theta_obser_arr),[length(phi_obser)],[1])]; % Monostatic observation angle
    theta = acos(r_hat(:,3));
    phi = mod(atan2(r_hat(:,2),r_hat(:,1)),2*pi);
    %phi = phi*(phi>=0) + (2*pi+phi)*(phi<0);
    sig_dir = sin(theta_r)*cos(phi_r)*(r*sin(theta).*cos(phi)-QS_arr(:,1)) + sin(theta_r)*sin(phi_r)*(r*sin(theta).*sin(phi)-QS_arr(:,2))...
    + cos(theta_r)*(r*cos(theta)-QS_arr(:,3)); 
    inv_Gamma_r = 1/(Gamma_r(1)*Gamma_r(4)-Gamma_r(2)*Gamma_r(3))*[Gamma_r(4),-Gamma_r(2),-Gamma_r(3),Gamma_r(1)]; %Inverse of Gamma_r
    omega = [sin(theta_obser_arr)*cos(theta_r)*cos(phi_obser-phi_r)-cos(theta_obser_arr)*sin(theta_r),...
         sin(theta_obser_arr)*sin(phi_obser-phi_r)]; 
    %omega = [repelem((theta_obser_arr - theta_r),[length(phi_obser)],[1]), sin(theta_obser_arr)*(phi_obser-phi_r)];
    ff_GB_decay = exp(1j*k/2*(omega(:,1).^2*inv_Gamma_r(1) + 2*omega(:,1).*omega(:,2)*inv_Gamma_r(2)+ ...
        omega(:,2).^2*inv_Gamma_r(4)));% Far field Gaussian Beam angular decay term 
    ff_phase_shift = exp(1j*k*dot(r_hat,QS_arr,2)); % Far field fraunhofer phase shift term
    Psi_thet_phi_mnst = A_Q*exp(-1j*k*l_b)*sqrt(1/((Gamma_r(1)*Gamma_r(4)-Gamma_r(2)*Gamma_r(3))))*ff_GB_decay.*ff_phase_shift.*(sig_dir>=del*r); %only considering the beams that 
    %produces forward contribution
    u_bistat_obser_lit  = u_bistat_obser_lit  + a_mu*transpose(Psi_thet_phi_mnst);
end
toc;
%% RCS Plotting 
Fig1 = figure(1);
plot((theta_obser_arr*180/pi),ones(1,length(theta_obser_arr)),'b-.','linewidth',2);
hold on
plot(theta_obser_arr*180/pi,4/a^2*abs(u_bistat_obser_lit).^2,'r','linewidth',2);
hold on 
ylabel('$\frac{\sigma}{\pi a^2}$','Interpreter','latex','FontSize',30);
xlabel('\theta [deg]');
set(gca,'FontSize',18); 
xlim([45,180]);
ylim([0,1.2]);
grid on
xlabel('\theta [deg]');
set(gca,'FontSize',18);
ylabel('$\frac{\sigma}{\pi a^2}$','Interpreter','latex','FontSize',30);
%saveas(Fig1,'Sphsim-SphRCS-a01b01k-x0z0a-nu0.20-bistat.png');
%% Error Plotting 
Fig2 = figure(2);
plot(theta_obser_arr*180/pi,20*log10(abs(ones(1,length(theta_obser_arr))-4/a^2*abs(u_bistat_obser_lit).^2)),'k','linewidth',2);
xlim([45,180]);
ylim([-100,0]);
grid on
xlabel('\theta [deg]');
set(gca,'FontSize',18);
ylabel('$\mathcal{E}$','Interpreter','latex');
saveas(Fig2,'Sphsim-Sphscatt-a03b1.5k-x0z0a-nu0.20-bistat.png');
