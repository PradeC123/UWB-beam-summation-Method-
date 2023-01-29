%%Date - 18/07/2021 
%% Far zone analysis
r = 10^100;
k = 1; b = 500; v_max = 0.33; % Expansion Parameter
n_1 = 5; n_2 = 5;
xi_bar = sqrt(2*pi*v_max/(k*b));
theta_b = acos(1-(n_1*xi_bar)^2-(n_2*xi_bar)^2);
phi_b = mod(atan2(n_2,n_1),2*pi);
sigma_i = [sin(theta_b)*cos(phi_b),sin(theta_b)*sin(phi_b),cos(theta_b)];
nu = [0,0,-1];
sigma_r = sigma_i ;%- 2*(dot(sigma_i,nu))*nu;
theta_r = acos(sigma_r(3));
phi_r = mod(atan2(sigma_r(2),sigma_r(1)),2*pi);
%% Far zone analysis (The field in the far zone) 
theta = 0:0.0005:pi;% Evevation angle.
%phi = 0:0.0005:2*pi;
%phi = pi;
phi = pi;
%[Theta,Phi] = meshgrid(theta,phi); 
%phi = 0;
sigma_r = sin(theta).*sin(theta_r).*cos(phi-phi_r)+cos(theta)*cos(theta_r);
eta_r1 = sin(theta)*cos(theta_r)*cos(phi-phi_r)-cos(theta)*sin(theta_r);
eta_r2 = sin(theta).*sin(phi-phi_r);
%% Far zone analysis (Simplified Expression)
eta_r11 = theta-theta_r;
eta_r22 = (repelem(sin(theta_r)*(phi-phi_r),length(theta)));
%% Far zone analysis
%theta = 0:0.001*pi/360:pi;% Evevation angle.
%phi = 0:0.01*2*pi/360:2*pi; 
rho_i = 1000;
Gamma_inc = [1/(rho_i+1j*b*cos(theta_b)^2),0;0,1/(rho_i+1j*b)];
theta_i = 0; %Angle of incidences on the sphere
a = 1000; %Radius of the curvature 
Gamma_Q = Gamma_inc + [2*cos(theta_i)/a,0;0,2/(a*cos(theta_i))];
inv_Gamma_r = 1/(Gamma_Q(1)*Gamma_Q(4)-Gamma_Q(2)*Gamma_Q(3))*[Gamma_Q(4),-Gamma_Q(2),-Gamma_Q(3),Gamma_Q(1)]; 
%% Field analysis
inv_Gamma_r = 1/(Gamma_Q(1)*Gamma_Q(4)-Gamma_Q(2)*Gamma_Q(3))*[Gamma_Q(4),-Gamma_Q(2),-Gamma_Q(3),Gamma_Q(1)]; 
F_1 = exp(1j*k/2*(inv_Gamma_r(1)*eta_r1(:).^2 + 2*inv_Gamma_r(2)*eta_r1(:).*eta_r2(:) +inv_Gamma_r(4)*eta_r2(:).^2)).*transpose((sigma_r>=1/2));
F_2 = exp(1j*k/2*(inv_Gamma_r(1)*eta_r11(:).^2+ 2*inv_Gamma_r(2)*eta_r1(:).*eta_r22(:)+ inv_Gamma_r(4)*eta_r22(:).^2));
%% Plotting the fields 
plot(theta*180/pi,20*log10(abs(F_2)),'r-.');
hold on 
plot(theta*180/pi,20*log10(abs(F_1)),'b');
ylim([-80 0])
%% 
QS_arr = [0,0,0];
theta = 0:0.0005:pi;% Evevation angle.
phi_array = 0:0.001:2*pi;
u_bistatic_1 = zeros(length(phi_array),length(theta));
u_bistatic_2 = zeros(length(phi_array),length(theta));
for cnt_i = 1:1:length(phi_array)
    disp(cnt_i);
    phi = phi_array(cnt_i);
    sig_dir = sin(theta_r)*cos(phi_r)*(r*sin(theta).*cos(phi)-QS_arr(1,1)) + sin(theta_r)*sin(phi_r)*(r*sin(theta).*sin(phi)-QS_arr(1,2))...
    + cos(theta_r)*(r*cos(theta)-QS_arr(1,3));
    eta_r1 = sin(theta)*cos(theta_r)*cos(phi-phi_r)-cos(theta)*sin(theta_r);
    eta_r2 = sin(theta).*sin(phi-phi_r);
    %eta_r11 = (theta-theta_r)-sin(theta)*cos(theta_r)*(phi-phi_r)^2/2;
    %eta_r22 = sin(theta_r)*(phi-phi_r)+cos(theta_r)*(phi-phi_r)*(theta-theta_r);
    eta_r11 = theta-theta_r; 
    eta_r22 = (repelem(sin(theta_r)*(phi-phi_r),length(theta)));
    %rho_i = 2000;
    %Gamma_Q = [1/(rho_i+1j*b*cos(theta_b)^2),0;0,1/(rho_i+1j*b)];
    %theta_i = pi/4; %Angle of incidences on the sphere
    %a = 2000; %Radius of the curvature 
    %Gamma_Q = Gamma_Q + [2*cos(theta_i)/a,0;0,2/(a*cos(theta_i))];
    %inv_Gamma_r = 1/(Gamma_Q(1)*Gamma_Q(4)-Gamma_Q(2)*Gamma_Q(3))*[Gamma_Q(4),-Gamma_Q(2),-Gamma_Q(3),Gamma_Q(1)];
    F_1 = exp(-1*k/2*imag(inv_Gamma_r(1)*eta_r1(:).^2 + 2*inv_Gamma_r(2)*eta_r1(:).*eta_r2(:) +inv_Gamma_r(4)*eta_r2(:).^2)).*transpose((sig_dir>=0*r));
    F_2 = exp(-1*k/2*imag(inv_Gamma_r(1)*eta_r11(:).^2 + 2*inv_Gamma_r(2)*eta_r11(:).*eta_r22(:)+ inv_Gamma_r(4)*eta_r22(:).^2));
    u_bistatic_1(cnt_i,:) = (F_1);
    u_bistatic_2(cnt_i,:) = (F_2);
end
plot(0)
hold on
imagesc(theta*180/pi,phi_array*180/pi,abs(u_bistatic_2-u_bistatic_1));
colorbar;
xlabel('\theta[deg]');
ylabel('\phi[deg]');
xlim([0,180]);
ylim([0,360])