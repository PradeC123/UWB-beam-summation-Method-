%% Chapter 4 PO numerical integration for sphere 
tic;
k = 1;
a = 1000;
del_theta = 2*pi/(10*k*a);
del_phi = 2*pi/(4*k*a);%% Bistatic RCS 
N_theta = 1000;
N_phi = 4000;
k = 1;
theta_arr = pi/2+(1/2)/N_theta:del_theta:pi/2*(2*N_theta-1/2)/N_theta;
phi_arr = 0:del_phi:2*pi;
theta_obs = pi;
phi_obs = 0;
[Theta_arr,Phi_arr] = meshgrid(theta_arr,phi_arr);
F = sin(Theta_arr).*(sin(Theta_arr).*sin(theta_obs).*cos(Phi_arr-phi_obs)...
    +cos(Theta_arr).*cos(theta_obs))...
    .*exp(-1j*k*a*((-sin(Theta_arr).*sin(theta_obs)...
    .*cos(Phi_arr-phi_obs)-cos(Theta_arr).*cos(theta_obs))+cos(Theta_arr)));
I = trapz(phi_arr,trapz(theta_arr,F,2));
Sph_mnst = (2*a^2*1j/(4*pi)*I);
Sph_mnst_RCS = 4*a^2*abs((2*1j/(4*pi)*I))^2;
toc;
%%
a = 1000;
N_theta = 4000;
N_phi = 16000;
k = 1;
del_theta = 2*pi/(10*k*a);
del_phi = 2*pi/(4*k*a);%% Bistatic RCS 
theta_obs_arr = 0:0.0001:180/180*pi;
N_theta = 3000;
N_phi = 9000;
k = 1;
del_theta = 2*pi/(10*k*a);
del_phi = 2*pi/(4*k*a);
theta_arr = pi/2+(1/2)/N_theta:del_theta:pi/2*(2*N_theta-1/2)/N_theta;
phi_arr = 0:del_phi:2*pi;
[Theta_arr,Phi_arr] = meshgrid(theta_arr,phi_arr);
U_bistat_1 = zeros(1,length(theta_obs_arr));
phi_obs = 0;
for cnt1 = 1:1:length(theta_obs_arr)
    theta_obs = theta_obs_arr(cnt1);
    F = sin(Theta_arr).*(sin(Theta_arr).*sin(theta_obs).*cos(Phi_arr-phi_obs)...
        +cos(Theta_arr).*cos(theta_obs))...
        .*exp(-1j*k*a*(-(sin(Theta_arr).*sin(theta_obs)...
    .*cos(Phi_arr-phi_obs)+cos(Theta_arr).*cos(theta_obs))+cos(Theta_arr)));
    I = trapz(phi_arr,trapz(theta_arr,F,2));
    Sph_mnst = (2*a^2*1j/(4*pi)*I);
    U_bistat_1(cnt1) = Sph_mnst;
end
toc;
% 

%%
tic;
fun = @(x,y) exp(-x.^2-y.^2);
ymax = @(x) sqrt(a^2 - x.^2);
ymin = @(x) -sqrt(a^2 - x.^2);
Q = quad2d(fun,-a,a,ymin,ymax);
toc;
tic;
fun = @(x,y) exp(-x.^2-y.^2);
ymax =  100*a;
ymin = -100*a;
Q = quad2d(fun,-100*a,100*a,ymin,ymax);
toc;
n_1 = 0;
n_2 = 0;
m_1 = 28;
m_2 = 0;
k = 100;
k_max = 100;
a = 10;
b = 10;
nu_max = 0.20;
x_bar = sqrt(2*pi*nu_max*b/k);
xi_bar = sqrt(2*pi*nu_max/(k*b));
theta_i = 0;
phi_i = 0;
tic;
conc = @(x,y) exp(-k*(x-m_1*x_bar).^2/(2*b) - k*(y-m_2*x_bar).^2/(2*b)).*exp(1j*k*(n_1*xi_bar-sin(theta_i)*cos(phi_i))*x)...
    .*exp(1j*k*(n_2*xi_bar-sin(theta_i)*sin(phi_i))*y).*exp(-1j*k*(m_1*x_bar*n_1*xi_bar+m_2*x_bar*n_2*xi_bar));
ymax = @(x) sqrt(a^2 - x.^2);
ymin = @(x) -sqrt(a^2 - x.^2);
nu_max^2*k^3/(pi*b*k_max^2)*integral2(conc,-a,a,ymin,ymax);
toc;

%%

clc; clear all ;
M = 1000 ;
N = 1000 ;
R1 = 0; % inner radius 
R2 = 10 ;  % outer radius
nR = linspace(R1,R2,M) ;
nT = linspace(0,2*pi,N) ;
%nT = pi/180*(0:NT:theta) ;
[R, T] = meshgrid(nR,nT) ;
% Convert grid to cartesian coordintes
X = R.*cos(T); 
Y = R.*sin(T);
[m,n]=size(X);
% Plot grid
figure
set(gcf,'color','w') ;
axis equal
axis off
box on
hold on
% Plot internal grid lines
for i=1:m
    plot(X(i,:),Y(i,:),'k','linewidth',1.5); 
end
for j=1:n
    plot(X(:,j),Y(:,j),'k','linewidth',1.5); 
end
%% Diffraction expansion coefficients 
n_1 = 0;
n_2 = 0;
m_1 = 30;
m_2 = 0;
k = 300;
k_max = 300;
a = 10;
b = 10;
nu_max = 0.20;
x_bar = sqrt(2*pi*nu_max*b/k);
xi_bar = sqrt(2*pi*nu_max/(k*b));
theta_i = 0;
phi_i = 0;
tic;
r_arr = linspace(0,a,a*20);
theta_arr = linspace(0,2*pi,a*20);
[R_arr,Theta_arr] = meshgrid(r_arr,theta_arr);
F = nu_max^2*k^3/(pi*b*k_max^2)*R_arr.*exp(-k/(2*b)*(R_arr.*cos(Theta_arr)-m_1*x_bar).^2)...
    .*exp(-k/(2*b)*(R_arr.*sin(Theta_arr)-m_2*x_bar).^2).*exp(1j*k*(n_1*xi_bar-sin(theta_i).*cos(phi_i))*R_arr.*cos(Theta_arr))...
    .*exp(1j*k*(n_2*xi_bar-sin(theta_i).*sin(phi_i))*R_arr.*sin(Theta_arr)).*exp(-1j*k*m_1*x_bar*n_1*xi_bar)...
    .*exp(-1j*k*m_2*x_bar*n_2*xi_bar);
trapz(theta_arr,trapz(r_arr,F));
toc;
%% Diffraction expansion coefficients Parameters that shall be computed 
k = 100;
k_max = 100;
a = 10;
b = 10;
nu_max = 0.20;
nu = 0.20;
x_bar = sqrt(2*pi*nu_max*b/k);
xi_bar = sqrt(2*pi*nu_max/(k*b));
theta_i = 0;
phi_i = 0;
M_1_diff = floor((a+2*x_bar)/x_bar);
m_1_diff = -M_1_diff:1:M_1_diff;
m_2_diff = -M_1_diff:1:M_1_diff;
%n_1max = round(1/xi_bar*sin(theta_inc)*cos(phi_inc));%
%n_2max = round(1/xi_bar*sin(theta_inc)*sin(phi_inc));
%C = 64; %Optimization parameter
%del_n = round(sqrt(C*k_max/(pi*nu*k)));
%n_1 = n_1max - del_n:1:n_1max + del_n;
%n_2 = n_2max -del_n:1:n_2max + del_n;
N_diff = floor(sin(30*pi/180)/xi_bar);
n_1_diff = -N_diff:1:N_diff;
n_2_diff = n_1_diff;
[N_1,N_2,M_1,M_2] = ndgrid(n_1_diff,n_2_diff,m_1_diff,m_2_diff);
Z_diff = [N_1(:),N_2(:),M_1(:),M_2(:)];
Z_diff(find((Z_diff(:,1)*xi_bar).^2 + (Z_diff(:,2)*xi_bar).^2 >1),:) = [];% Elimination of the beams propagating in evanescent mode
Z_diff(find(sqrt((Z_diff(:,3)*x_bar).^2 + (Z_diff(:,4)*x_bar).^2)>(a+2*x_bar)),:)=[];
Z_diff(find(sqrt((Z_diff(:,3)*x_bar).^2 + (Z_diff(:,4)*x_bar).^2)<(a-2*x_bar)),:)=[];
r_arr = linspace(0,a,a*30); %Radial parameter for integration 
theta_arr = linspace(0,2*pi,a*30);%Azimuthal parameter for integration 
[R_arr,Theta_arr] = meshgrid(r_arr,theta_arr);
for cnt1 = 1:1:length(Z_diff)
    n_1 = Z_diff(cnt1,1);
    n_2 = Z_diff(cnt1,2);
    m_1 = Z_diff(cnt1,3);
    m_2 = Z_diff(cnt1,4);
    F = nu_max^2*k^3/(pi*b*k_max^2)*R_arr.*exp(-k/(2*b)*(R_arr.*cos(Theta_arr)-m_1*x_bar).^2)...
    .*exp(-k/(2*b)*(R_arr.*sin(Theta_arr)-m_2*x_bar).^2).*exp(1j*k*(n_1*xi_bar-sin(theta_i).*cos(phi_i))*R_arr.*cos(Theta_arr))...
    .*exp(1j*k*(n_2*xi_bar-sin(theta_i).*sin(phi_i))*R_arr.*sin(Theta_arr)).*exp(-1j*k*m_1*x_bar*n_1*xi_bar)...
    .*exp(-1j*k*m_2*x_bar*n_2*xi_bar);
    Z_diff(cnt1,5) = trapz(theta_arr,trapz(r_arr,F));
end

