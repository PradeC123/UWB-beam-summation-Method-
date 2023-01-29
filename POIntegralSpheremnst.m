%% PO Code
function [u_scatt_theta, t] = POIntegralSpheremnst(k_min,k_max,del_k,K)
    a = 1;
    k_arr = linspace(k_min,k_max,del_k);
    ka = max(k_arr);
    %ka_max = max(k_arr);
    k = ka; %Our Test case
    theta_obs_ang = 180/180*pi; %PRANAAAAV
    %theta_obs_ang = pi;
    u_scatt_theta = zeros(1,length(theta_obs_ang));
    u_theta_summsnt = [];
    theta_m_cnst = deg2rad(179); %angle after which the phi distribution on the shell are constant 
    N_theta = floor(K*ka);
    del_thet = (pi/2)/N_theta;
    phi_obs = pi;
    theta_0 = sqrt(pi)/sqrt(min(k_arr)*a);
    filt_1 = [];
    theta_m1 = [];
    for cnt1 = 1:1:length(theta_obs_ang)
        tic;
        theta_obs_angle_a = theta_obs_ang(cnt1);
        sum_bist = 0;
        for cnt2 = 1:1:N_theta
            theta_m = (N_theta+1-cnt2-1/2)*del_thet + pi/2;
            if theta_m >= theta_m_cnst 
                N_phi_m = floor(2*K*ka*sin(theta_m_cnst));
            else
                N_phi_m = floor(2*K*ka*sin(theta_m));
            end
            if (theta_m >= pi/2)&&(theta_m <= pi/2 + 5*pi/180) %Lower end point filter (Tapper_L) 
                tapper_m = -2*(180/(5*pi))^3*(theta_m-pi/2)^3 + 3*(180/(5*pi))^2*(theta_m-pi/2)^2;
            %elseif (theta_m >= pi - 5*pi/180)&&(theta_m <= pi) 
             %   tapper_m = 2*(180/(5*pi))^3*(theta_m-pi)^3 + 3*(180/(5*pi))^2*(theta_m-pi)^2;
            else
                tapper_m = 1; 
            end
            k_theta = 2;
            k_min = min(k_arr);
            sigma_1 = sqrt(2*pi)/(sqrt(k_min*a));
            if (theta_m>pi-k_theta*theta_0)
                filt = 1;
            else
                filt = exp(-(theta_m-(pi-k_theta*theta_0)).^2/(2*sigma_1^2));
            end
            %filt = 1;
            filt_1 = [filt_1  filt];
            %K_tild = 0.01;
            %K_theta = 400;
            %c = -(ka)*log(K_tild)/(K_theta)^2;
            %sigma_1 = sqrt(K_theta/10)*sqrt(pi)/(sqrt(ka));
            %end
            %filt = exp(-(theta_m-pi).^2/(2*sigma_1^2)); % Fitler to filter out the extra oscillations
            del_phi_m = 2*pi/N_phi_m;
            phi_nm_arr = linspace(0+rand(1),2*pi+rand(1),N_phi_m);
            f_thetn_phim = (2*1j*a^2)/(4*pi)*(sin(theta_m)*sin(theta_obs_angle_a)*cos(phi_obs-phi_nm_arr)+...
                    cos(theta_obs_angle_a)*cos(theta_m))*sin(theta_m)*tapper_m;
            Phi_thetm = -a*sin(theta_m)*sin(theta_obs_angle_a)*cos(phi_obs-phi_nm_arr)...
                -a*cos(theta_obs_angle_a)*cos(theta_m) + a*cos(theta_m);
            %S_m_theta = k*(f_thetn_phim.*exp(-1j*k*Phi_thetm)*del_phi_m*del_thet)*filt;
            S_m_theta = k_arr.*(transpose(f_thetn_phim).*exp(-1j*k_arr.*transpose(Phi_thetm)).*del_phi_m*del_thet)*filt;
            sum_bist = sum_bist + sum(S_m_theta);
            %u_theta_summsnt = [u_theta_summsnt;4/a^2*abs(sum_bist).^2];
            %u_theta_summsnt = [u_theta_summsnt,abs(sum_bist)];
            %u_scatt_theta(cnt1) = sum_bist;
            theta_m1 = [theta_m1 theta_m];
            %if theta_m < pi-8*theta_0
            %    break 
            %end
            if filt <= 10^-6
                break 
            end
        end
        t = toc;
    end
end
%% 



