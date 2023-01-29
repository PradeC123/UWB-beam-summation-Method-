  %% Chapter 4 Part 4: RCS at an particular observation point
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
function [RCS_obser,beamseps1,beamseps2] = Chapter4RCSobser(Z,theta_obser,phi_obser,a,a_mu,k,eps1,eps2)
    Sph_mnst = 0;
    Z_r = Z;
    Z_r(:,13) = a_mu;
    r = 10000000000000*a;%Far Zone
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
        if (sig_dir>=0.85*r) %In the lit region
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
    RCS_obser = (Sph_mnst);
    beamseps1 = sum(abs(Z_r(:,14))>eps1);
    RCS_obser_eps1 = 4/(a^2)*abs(sum(Z_r(abs(Z_r(:,14))>eps1,14)))^2;
    beamseps2 = sum(abs(Z_r(:,14))>eps2);
    RCS_obser_eps2 = 4/(a^2)*abs(sum(Z_r(abs(Z_r(:,14))>eps2,14)))^2;
    %disp(RCS_obser_eps1);
    %disp(RCS_obser_eps2);
end




