%% Computing the number of significant GBs needed to sum up obtain monostatic RCS
k_arr =  linspace(250,8000,100);%frequency band
no_of_beams = zeros(1,length(k_arr)); %No of beams needed to be summed
RCS_tracker = zeros(1,length(k_arr)); %Tracks the RCS for a particular frequency
for cntk = 1:1:length(k_arr) 
    theta_inc = 0; %inicident angle of propagation of the field
    phi_inc = 0; %incident angle of propagation of the field
    k_max = 1; % The normalized wavenumber also the highest frequency in the frequency band
    k = 1;  %k = 0.5 0.75 1 observation frequency
    a = k_arr(cntk);% The radius of the sphere
    b = 4/4*a; %b = 1000,2000,4000, Collimation distance in the frequency normalized system
    W = sqrt(b/k);
    nu = 0.20;
    x_bar = sqrt(2*pi*nu*b/k_max);%Frequency independent lattice
    xi_bar = sqrt(2*pi*nu/(b*k_max));
    Theta_d1 = sqrt(1/(k*b));
    z_b = -a; %location of the expansion plane
    %L = x1_max; %the width of the square aperture in x1 axis. Onl relevant for truncation of the beams
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Respective indices in m_1,m_2,n_1,n_2
    M1 = floor(1.05*a/x_bar);
    M2 = floor(1.05*a/x_bar);
    m_1 = -M1:1:M1;
    m_2 = m_1;%
    n_1max = round(1/xi_bar*sin(theta_inc)*cos(phi_inc));%
    n_2max = round(1/xi_bar*sin(theta_inc)*sin(phi_inc));
    C = 16; %Optimization parameter
    del_n = round(sqrt(C*k_max/(pi*nu*k)));
    n_1 = n_1max - del_n:1:n_1max + del_n;
    n_2 = n_2max -del_n:1:n_2max + del_n;
    %N = floor(sin(30*pi/180)/xi_bar);
    %n_1 = -N:1:N;
    %n_2 = n_1;
    [N_1,N_2,M_1,M_2] = ndgrid(n_1,n_2,m_1,m_2);
    Z = [N_1(:),N_2(:),M_1(:),M_2(:)];% Generates all the possible combination of n1,n2,m1,m2
    Z(find((Z(:,1)*xi_bar).^2 + (Z(:,2)*xi_bar).^2 >1),:) = [];% Elimination of the beams propagating in evanescent mode
    theta_n = acos(sqrt(1-xi_bar^2*(Z(:,2).^2+Z(:,1).^2)));
    phi_n = atan2(Z(:,2),Z(:,1));
    phi_n = phi_n .* (phi_n >= 0) + (phi_n + 2 * pi) .* (phi_n < 0);
    Z(:,5) = z_b;
    Z(:,6) = theta_n;
    Z(:,7) = phi_n;
    m1 = Z(:,3);
    m2 = Z(:,4);
    n1 = Z(:,1);
    n2 = Z(:,2);
    % Neglecting the edge effects[difraction manifolds]...
    a_mu = nu^2*2*(k/k_max)^2*exp(-k*b*(sin(theta_inc)*cos(phi_inc)-Z(:,1)*xi_bar).^2/2).*exp(-k*b*(sin(theta_inc)*sin(phi_inc)-Z(:,2)*xi_bar).^2/2).*...
        exp(-1j*k*sin(theta_inc)*cos(phi_inc)*Z(:,3)*x_bar).*exp(-1j*k*sin(theta_inc)*sin(phi_inc)*Z(:,4)*x_bar).*exp(-1j*k*z_b*cos(theta_inc));
    % a_mu = Rec_coeff_PW(theta_inc,phi_inc,nu,b,k,m1,m2,n1,n2,z_b);
    Z(:,8) = a_mu;
    %Z(find((abs(a_mu))<10^-15),:)=[];
    Z(find(abs(a_mu)<10^-6),:)=[];
    Z_r = zeros(length(Z),11);
    QS = zeros(length(Z),3);
    theta_r = zeros(length(Z),1);
    phi_r = zeros(length(Z),1);
    l_b = zeros(length(Z),1);
    hit = zeros(length(Z),1);
    Gamma_r = zeros(length(Z),4);
    A_Q = zeros(length(Z),1);
    c1 = [0,0,0];
    r1 = a;
    %c2 = [-3*a,0,0];
    %r2 = a;
    % [X1,Y1,Z1] = sphere;
    % r = 200;
    % X2 = X1 * r;
    % Y2 = Y1 * r;
    % Z2 = Z1 * r;
    % surf(X2,Y2,Z2);
    % hold on;
    for cnt1 = 1:1:length(Z)
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
    %Z_r(find(Z_r(:,13)<10^-10),:)=[];
    %Z_r(find((abs(Z_r(:,13))/max(abs(Z_r(:,13)))10^-10)),:)=[];
    Sph_mnst = 0;
    m = (1/2*k*a)^(1/2);
    theta_obser = 180*pi/180 - theta_inc;
    phi_obser = 0*pi + phi_inc;
    r = 100*a;
    for cnt1 = 1:1:length(Z_r)
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
        %phi = phi*(phi>=0) + (2*pi+phi)*(phi<0);
        %sig_dir = sin(theta_r)*cos(phi_r)*sin(theta)*cos(phi) + sin(theta_r)*sin(phi_r)*sin(theta)*sin(phi)...
            %+ cos(theta_r)*(cos(theta));
        sig_dir = sin(theta_r)*cos(phi_r)*(r*sin(theta)*cos(phi)-QS(1)) + sin(theta_r)*sin(phi_r)*(r*sin(theta)*sin(phi)-QS(2))...
            + cos(theta_r)*(r*cos(theta)-QS(3));
        if (sig_dir>=0) %In the lit region
            omega = [sin(theta_obser)*cos(theta_r)*cos(phi_obser-phi_r)-cos(theta_obser)*sin(theta_r);sin(theta_obser)*sin(phi_obser-phi_r)];
            %omega = [(theta-theta_r); sin(theta_r)*(phi-phi_r)];
            %omega = [sin(theta_obser)*cos(phi_obser)- sin(theta_r)*cos(phi_r);sin(theta_obser)*sin(phi_obser)- sin(theta_r)*sin(phi_r)];
            ff_GB_decay = exp(1j*k/2*transpose(omega)*((Gamma_r)^-1)*omega);% Far field Gaussian Beam angular decay term
            ff_phase_shift = exp(1j*k*dot(r_hat,QS)); % Far field fraunhofer phase shift term
            Psi_thet_phi_mnst = A_Q*sqrt(1/(det(Gamma_r)))*ff_GB_decay*ff_phase_shift;
            Sph_mnst = Sph_mnst + a_mu*Psi_thet_phi_mnst*exp(-1j*k*l_b);
            Z_r(cnt1,14) = a_mu*Psi_thet_phi_mnst*exp(-1j*k*l_b);
        end
    end
    no_of_beams(cntk) = length(find(abs(Z_r(:,14))>10^-3));
    RCS_tracker(cntk) = 4/a^2*abs(sum(Z_r(find(abs(Z_r(:,14))>10^-3),14))).^2;
end
%% Plotting the results 
semilogy(log2(f_arr),(k_arr).^2,'k','linewidth',2);
hold on;
semilogy(log2(f_arr),no_of_beams,'r','linewidth',2);
xticks([2,3,4,5,6,7]);
xlabel('f(Ghz)','FontSize',15);
ylabel("N");
set(gca, 'XTickLabel',[])                      %# suppress current x-labels
xt = get(gca, 'XTick');
yl = get(gca, 'YLim');
str = cellstr( num2str(xt(:),'2^{%d}') );      %# format x-ticks as 2^{xx}
hTxt = text(xt, yl(ones(size(xt))), str, ...   %# create text at same locations
'Interpreter','tex', ...                   %# specify tex interpreter
'VerticalAlignment','top', ...             %# v-align to be underneath
'HorizontalAlignment','center');
ylim([10^3,10^9]);
legend('PO Integral','GBSM');