% Date - 02/09/2020 (18:42) Israel Daylight Time (Started on)
tic;
%Initial beam properties.
k_max = 1;
v_max = 0.2546;
k = 1; % the wavenumber of the inital beam.
%Target Properties(Kept in the free space)
a = 500; %radius of the curvature. 
%The optimized way to find the beam width such that the phase error between
%the incident and the reflected beam is less than pi/4
%r = [-8*k,4*k*a,pi*a,-pi*a^2/4]; 
%q = 2; %q is the tunning parameter to decide for the phase error measurement q beamwidths away from te beam axis.
%W = sqrt(min(abs(roots(r)))*2/(q^2*a))*a;% the beam width of the beam.
%W = 20*lambda;
% -----------------------------------------------------------------------
%b = 100;% The Rayleighs length of the inital beam.
%b = 40;
%W = sqrt(b/k);
%b = 1;
b =  500;
W = sqrt(b/k);
theta_inc = 0;
phi_inc = 0;
Thet_d1 = 1/sqrt(k*b*cos(theta_inc)^2);
Thet_d2 = 1/sqrt(k*b);
% Choosing the frame parameters
nu_max = 0.254;
nu = nu_max*k/k_max;
x_bar = sqrt(2*pi*nu_max*b/k_max);
xi_bar = sqrt(2*pi*nu_max/(b*k_max));
z_b = 0;
%M = ceil((2*a/(x_bar))); %discretization a/x_bar in addition to that 5 beamwidth addition.
M = ceil((2*a/(x_bar)));
N = 20;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
m_1 = -M:1:M;
m_2 = 0;
n_1 = -N:1:N;
n_2 = 0;
z_Q = z_b;% The phase shift of the beam.
[N_1,N_2,M_1,M_2] = ndgrid(n_1,n_2,m_1,m_2);
Z = [N_1(:),N_2(:),M_1(:),M_2(:)];% Generates all the possible combination of n1,n2,m1,m2,zQ
Z(find((Z(:,1)*xi_bar).^2 + (Z(:,2)*xi_bar).^2 >1),:) = []; %Discarding beams in evanescent mode.
theta_n = acos(sqrt(1-xi_bar^2*(Z(:,2).^2+Z(:,1).^2)));
phi_n = atan2(Z(:,2),Z(:,1));
phi_n = phi_n.* (phi_n >= 0) + (phi_n + 2 * pi) .* (phi_n < 0);
Z(:,5) = z_b; 
Z(:,6) = theta_n;
Z(:,7) = phi_n; 
m1 = Z(:,3);
m2 = Z(:,4);
n1 = Z(:,1);
n2 = Z(:,2);
a_mu = Rec_coeff_PW(theta_inc,phi_inc,nu,b,k,m1,m2,n1,n2,z_b);
Z(:,8) = a_mu;
%Z(find((abs(Z(:,8))/max(abs(Z(:,8)))<10^-6)),:)=[];
%Z(find((Z(:,3)==0).*(Z(:,4)==0).*(Z(:,2)==0).*(Z(:,1)==0)),:);% mu = 0
toc;
%% Computing the hit point on the sphere.
% center of the sphere 
c = [0,0,a];
c_Z = zeros(length(Z),3);
for i = 1:1:length(Z)
    c_Z(i,:) = c;
end
X_0 = [Z(:,3)*x_bar,Z(:,4)*x_bar,Z(:,5)];% all the relevant combination that contributes to the far zone.
sig_i = [sin(Z(:,6)).*cos(Z(:,7)),sin(Z(:,6)).*sin(Z(:,7)),cos(Z(:,6))];% the propagation axis of the beams
[QS,nu_vect,l_b,hit] = sph_intersection(X_0,sig_i,c_Z,a); 
%x_0 is the origin of the beam axis, QS is the point the lies on the sphere
%nu_vect is the outward normal to the point Q, l_b is the distance
%propagated by the GB, hit is the number of times a GB enounters a
%Scatterer.
QS(find(hit==0),:) = [];
nu_vect(find(hit==0),:) =[];
l_b(find(hit==0))=[];
X_0(find(hit==0),:)=[];
sig_i(find(hit==0),:)=[];
Z(find(hit==0),:)=[]; %Removing all the case that dont lie on the sphere.
hit(find(hit==0))=[];
%%
%Calculating the direction of the reflected field.
sig_r = (sig_i-(2*(dot(sig_i,nu_vect,2)).*nu_vect))./vecnorm(sig_i-2*(dot(sig_i,nu_vect,2)).*nu_vect,2,2);
theta_r = acos(sig_r(:,3));% Theta angle of reflected beam
phi_r = atan2(sig_r(:,2),sig_r(:,1)); % corresponding Phi angle of the reflected beam
phi_r = phi_r.* (phi_r >= 0) + (phi_r + 2 * pi) .* (phi_r < 0);
% pi-theta_inc-5*Thet_d1 monostatic angle.
%Z_0 = zeros(length(Z),3);
%%
Sph_mnst = 0;
for cnt1 = 1:1:length(Z_monstat(:,1))
    theta_GB = Z(cnt1,6);
    phi_GB = Z(cnt1,7);
    QS_1 = [QS(cnt1,1),QS(cnt1,2),QS(cnt1,3)];%This point lies on the spher
    Rot_etasig_xyz = [cos(theta_GB)*cos(phi_GB),cos(theta_GB)*sin(phi_GB),-sin(theta_GB);...
        -sin(phi_GB),cos(phi_GB),0;...
        sin(theta_GB)*cos(phi_GB),sin(theta_GB)*sin(phi_GB),cos(theta_GB)];
    nu_otnormal = nu_vect((cnt1),:);
    l_b_GB = l_b((cnt1),:);
    Gamma_GB = [1/(1j*b),0;0,1/(1j*b)];
    R_Q = QS_1;% in our case.
    [Gamma_new_Q_inc,Rot_etasig_new_inc,Rot_Q_xinu_xyz] = beam_rot_POI(Gamma_GB,Rot_etasig_xyz,nu_otnormal,l_b_GB);
    [Gamma_r,Rot_lbc_eta_sig_r,Amp_GB_r] = GB_Refl(Gamma_GB,Gamma_new_Q_inc,l_b_GB,Rot_etasig_new_inc,Rot_Q_xinu_xyz,a);
    [Psi_thet_phi_mnst,ff_GB_decay,ff_phase_shift] = GB_FF_monost(Gamma_r,Rot_lbc_eta_sig_r,Amp_GB_r,theta_inc,phi_inc,k,R_Q);
    Sph_mnst = Sph_mnst + Z(cnt1,8)*Psi_thet_phi_mnst*exp(-1j*k*l_b_GB);
    Z(cnt1,9) = Psi_thet_phi_mnst;
    Z(cnt1,10) = ff_phase_shift;
    %Z_0(cnt1,1) = Psi_thet_phi_mnst;
    %Z_0(cnt1,2) = ff_GB_decay;
    %Z_0(cnt1,3) =  ff_phase_shift;
end
disp(abs(Sph_mnst).^2*4*pi/(a^2*pi));

