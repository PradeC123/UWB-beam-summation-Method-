%% 3D Multi_intersection Calculator
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input of the function 

%k is the wavenumber of the operating frequency 

%x_0 is the origin of the beam 

%theta_i is the angle of elevation of the incident GB. 

%phi_i is the angle of azimuthal of the incident GB. 

%Gamma_i is the complex curvature matrix for the incident GB

%c1 is the center of the sphere 
%r1 is the radius of the sphere

%%One should add the cn and rn parameter according to number of spheres in
%%the target Domain.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output of the File
% QS is the exit point 
% theta_GB_r is the elevation angle of the reflected field exiting the
% target domain 

% phi_GB_r is the azimuthal angle ----'''------
% Gamma_r_Q is the Complex curvature matrix for the reflected field at the
% exit point 

% A_Q is the amplitude and the phase accumulation term for each GB
% interaction

%hit_tot is the number of times a GB interacted in the Target Domain
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [QS,theta_sig_r,phi_sig_r,l_b,Gamma_r_tot,A_Q,hit_tot] = multi_intersection_3D(x_0,theta_i,phi_i,Gamma_i,c1,r1)%,c2,r2)
%,c2,r2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 1: Incident Field charachteristics on the expansion plane;
sigma_i = [sin(theta_i)*cos(phi_i),sin(theta_i)*sin(phi_i),cos(theta_i)]; %direction of propagation of the incident GB
eta_i1 = [cos(theta_i)*cos(phi_i),cos(theta_i)*sin(phi_i),-sin(theta_i)]; %transvere direction of the GB 
eta_i2 = [-sin(phi_i),cos(phi_i),0]; 
Gamma_i = transpose(reshape(Gamma_i,[2,2])); %Complex curvature matrix of a GB 
tar_dmn = 1;% beam from the phase-space lattice are always in the target domain. 
hit = 0;% any hit on the surface shall trigger this parameter to one.
sum_lb = 0;% Keeping tracking of the distance propagates by the GB along its axis while interacting with multiple scatteres in the target domain
hit_tot = 0;% Computing the totat number of interaction of a GB in the target domain
A_i_Q = 1;% Initial amplitude
while (tar_dmn == 1)&&(hit_tot)<10 
    Sph_info = [c1,r1];%c2,r2];%;;c3,r3;c4,r4]; %First two elements are the center of the sphere and the last one in the radius of the sphere
    Delta1 = (dot(sigma_i,x_0-c1,2)).^2 - (vecnorm(x_0-c1,2,2).^2 -r1^2);
    %Delta2 = (dot(sigma_i,x_0-c2,2)).^2 - (vecnorm(x_0-c2,2,2).^2 -r2^2);
    %Delta3 = (dot(sigma_i,x_0-c3,2)).^2 - (vecnorm(x_0-c3,2,2).^2 -r3^2);
    %Delta4 = (dot(sigma_i,x_0-c4,2)).^2 - (vecnorm(x_0-c4,2,2).^2 -r4^2);
    Delta_sph = [Delta1];%,Delta2];%, Delta2]; %Delta3, Delta4]; %,Delta3, Delta4]; %WE must demand taht the delta must be+ for both intersections
    Delta_sph(Delta_sph>=0) = 1; %Only considering the 1/2 intersections with the spheres.
    Delta_sph(Delta_sph<-0.00000000000000000000001) = nan;% disregarding the solution with no intersection.
    vect_dir = (-1)*[ dot(sigma_i,x_0-c1,2)];%,dot(sigma_i,x_0-c2,2)];% dot(sigma_i,x_0-c3,2),dot(sigma_i,x_0-c4,2)];%We shall have intersections if any Positve 
    vect_dir(vect_dir<0) = nan; 
    intersect_sol = vect_dir.*Delta_sph;
    if isnan( min(intersect_sol)) == true 
        if hit == 0%no intersection with the sphere 
            c_r = [0,0,0];% Location of the target domain sphere
            r_domain = 10^5;% Location of the radius of the domain sphere
            Delta_1 = (dot(sigma_i,x_0-c_r,2)).^2 - (vecnorm(x_0-c_r,2,2).^2 -r_domain^2);
            d_1 = -dot(sigma_i,x_0-c_r,2) + sqrt(Delta_1);
            d_2 = -dot(sigma_i,x_0-c_r,2) - sqrt(Delta_1);
            if d_1 > d_2 
                d = d_2;
            else 
                d = d_1;  
            end
            tar_dmn = 0; %Emerges out of the target domain
            QS = x_0 + d.*sigma_i;
            nu_vect = (QS - r_domain)./vecnorm((QS - r_domain),2,2);
            sig_r_Q = sigma_i - 2*(dot(sigma_i,nu_vect))*nu_vect; %Computing the reflected field direction, the reflected field shall acts as the new incident GB for our algorithm
            sig_r = sig_r_Q/norm(sig_r_Q);
            theta_sig_r = theta_i; 
            phi_sig_r = phi_i;
            Gamma_r_tot = (reshape(transpose(Gamma_i),[1,4]));
            l_b = vecnorm(QS-x_0,2,2);    
            A_Q = 1;
            hit = 0; %no interaction of the GB from any sphere in the target domain
        else %The beam will now emerge in the scattering domain
%             d = 0:0.1:4000; 
%             X_r_path =  x_0 + transpose(d).*sigma_i; %Line equation for the reflected field axis.
%             plot(X_r_path(:,3),X_r_path(:,1),'k','linewidth',2 );
%             hold on;
%             %shg;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %Step 5: Rotation of the transverse coordinate 
            QS = QS_inter; %last hit before going to the farzone
            sigma_r = sigma_i; % propagating direction of the relfected beam in the scattering domain
            eta_r1 = eta_i1;% transvere eta_r1 direction of the relfected beam in the scattering domain
            eta_r2 = eta_i2; % transverse eta_r2 direction of the relfected beam in the scattering domain
            theta_sig_r = acos(sigma_r(3)); %reflected beam's elevation angle             
            phi_sig_r = mod(atan2(sigma_r(2),sigma_r(1)),2*pi); %reflected beams's azimuthal angle
            % Beam transverse coordinate rotation
            eta_tilde_1r = [cos(theta_sig_r)*cos(phi_sig_r),cos(theta_sig_r)*sin(phi_sig_r),-sin(theta_sig_r)];
            eta_tilde_2r = [-sin(phi_sig_r),cos(phi_sig_r),0];
            Phi_eta  = [dot(eta_tilde_1r,eta_r1),dot(eta_tilde_1r,eta_r2);dot(eta_tilde_2r,eta_r1),dot(eta_tilde_2r,eta_r2)];
            Gamma_tild_r = ((Phi_eta))*(Gamma_i)*(transpose(Phi_eta));
            Gamma_r_tot = (reshape(transpose(Gamma_tild_r),[1,4]));% Gamma_r at the exit point;
            A_Q =  A_i_Q; %Amplitude accumulated by GB at every interaction.
            l_b = sum_lb;% Distance propagated by the GB until the last intersection.
            hit = hit_tot; % Total number of hits the beam occured 
            tar_dmn = 0; %Emerges out of the target domain
        end
    else
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Step 2: Computation of the interaction of the GB on the boundary of the sphere.
        [~,true_intersect_indx] = min( intersect_sol); %computing the first interaction with the multiple scattering environment
        c_center = Sph_info(true_intersect_indx,1:3);% Center of the first intersection of the sphere 
        c_radius = Sph_info(true_intersect_indx,4); % Radius of the sphere of the first intersection of the sphere
        Delta_inter = (dot(sigma_i,x_0-c_center,2)).^2 - (vecnorm(x_0-c_center,2,2).^2 -c_radius^2);
        d_inter1 = -dot(sigma_i,x_0-c_center,2) + sqrt(Delta_inter);
        d_inter2 = -dot(sigma_i,x_0-c_center,2) - sqrt(Delta_inter);
        if d_inter1 > d_inter2 
            d_inter = d_inter2;% Minimum d_inter corresponds to the first interaction with the sphere
        else 
            d_inter = d_inter1; 
        end
        QS_inter = x_0 + d_inter.*sigma_i; %Specular point
        nu_vect_inter = (QS_inter - c_center)./vecnorm((QS_inter - c_center),2,2); %Normal to the plane of incidence 
        l_b_inter = vecnorm(QS_inter-x_0,2,2); %Distance from the origin to the specular point
        hit = 1;%We encounted the first hit
        tar_dmn = 1; %The GB is still inside the target domain
        sum_lb = sum_lb + l_b_inter; % sum of all the distance accumulated by GB from the origin
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Plotting the reflected beams in the space
%         Fig2 = figure(3);
%         [X,Y,Z] = sphere;
%         r = 1000;
%         X2 = X * r;
%         Y2 = Y * r;
%         Z2 = Z * r;
%         surf(X2,Y2,Z2);
%         hold on;
%         l_b_crr = 0:l_b_inter/1000:l_b_inter;
%         X_i_path = x_0+transpose(l_b_crr).*sigma_i; %Line equation of all the incident field.
%         plot3(X_i_path(:,1),X_i_path(:,2),X_i_path(:,3),'r','linewidth',2);
%         hold on;
%         d = 0:1:l_b_inter;
%         Nu_r_path = QS_inter + transpose(d).*nu_vect_inter;
%         plot3(Nu_r_path(:,1),Nu_r_path(:,2),Nu_r_path(:,3),'--k','linewidth',2);
%         hold on;
%         sig_r_Q = sigma_i - 2*(dot(sigma_i,nu_vect_inter))*nu_vect_inter; %Computing the reflected field direction, the reflected field shall acts as the new incident GB for our algorithm
%         sig_r_Q = sig_r_Q/norm(sig_r_Q);
%         X_r_path = QS_inter+transpose(l_b_crr).*sig_r_Q;
%         plot3(X_r_path(:,1),X_r_path(:,2),X_r_path(:,3),'b','linewidth',2);
%         shg;
%         pause(0.01);
%         hold on;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Sketching the incident vectors along w
        l_b_crr = 0:l_b_inter/1000:l_b_inter;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Step 3: Establishing a plane of incident and Rotation of the
        % incident GB in the plane of incidence
        Gamma_i_Q = inv(inv(Gamma_i)+l_b_inter*eye(2)); % Complex Curvatue matrix at the specular point
        A_i_Q = A_i_Q*sqrt(det(Gamma_i_Q)/det(Gamma_i));% The Complex Amplitude at the specular point
        xi_2 = cross(nu_vect_inter,-sigma_i,2); %the vector perpendicular to the plane of incidence
        if vecnorm(xi_2,2,2)==0 %This case implies Plane of incidence doesnot exists and nu_nect and sigma are colinear.
            xi_2 = eta_i2; 
            xi_1 = cross(xi_2,nu_vect_inter,2);
        else
            xi_2 = xi_2./vecnorm(xi_2,2,2);%normalization of the vector perpendicular to plane of incidence 
            xi_1 = cross(xi_2,nu_vect_inter,2); 
            xi_1 = xi_1./vecnorm(xi_1,2,2);% Normalizaton of the vector in the plane of incidence.
        end
        % Choosing the desired beam coordinate system such that, eta_new_i1 is in
        % the perpendicular direction of the plane of incidence and eta_new_i2 is in the plane of
        % incendence
        eta_new_i2 = xi_2;
        eta_new_i1 = cross(eta_new_i2,sigma_i);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % The matrix that rotates the plane of the transverse coordinate system to
        % the new desired plane such that it satisfies the specified conditions as
        % mentioned above.
        theta_eta_12 = [ dot(eta_new_i1,eta_i1,2), dot(eta_new_i1,eta_i2,2); ...
        dot(eta_new_i2,eta_i1,2), dot(eta_new_i2,eta_i2,2)];
        % Roation of the Gamma_i_Q in the plane of incidence
        Gamma_new_i_Q = theta_eta_12*Gamma_i_Q*transpose(theta_eta_12);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Step 4: Establishing a plane of incident and Rotation of the
        % incident GB in the plane of incidence
        theta_i = acos(-dot(sigma_i,nu_vect_inter)); % The initial angle of incidence
        % xi_1, xi_2 and nu are the local coordinate on the sphere, at the point Q
        % Establishing the matrix xi1xi2nu as a fucntion of the local
        % coordinate system
        Rot_lcs_xi_nu(1,:) = xi_1;
        Rot_lcs_xi_nu(2,:) = xi_2;
        Rot_lcs_xi_nu(3,:) = nu_vect_inter; 
        % Establishing the matrix (eta_new_i1,eta_new_i2,sigma_i) as a
        % fucntion of the global coordinate sytem 
        Rot_lbc_eta_sig_new_int(1,:) = eta_new_i1;
        Rot_lbc_eta_sig_new_int(2,:) = eta_new_i2;
        Rot_lbc_eta_sig_new_int(3,:) = sigma_i;
        % Matrix relation between (eta_i_1, eta_i_2, sigma_i) and (xi_1, xi_2, nu)
        Rot_etasig_nuxi = [dot(eta_new_i1,xi_1,2), dot(eta_new_i1,xi_2,2), dot(eta_new_i1,nu_vect_inter,2);...
        dot(eta_new_i2,xi_1,2),dot(eta_new_i2,xi_2,2),dot(eta_new_i2,nu_vect_inter,2);...
        dot(sigma_i,xi_1,2),dot(sigma_i,xi_2,2), dot(sigma_i,nu_vect_inter,2)];
        % Projection matrix of (eta_i_1,eta_i_2) on (xi_1,xi_2)
        Theta_eta_i = Rot_etasig_nuxi([1,2],[1,2]);
        % Matrix relation between (eta_r_1, eta_r_2, sigma_i) and (xi_1, xi_2, nu) 
        Rot_etsig_nu_r = [-dot(eta_new_i1,xi_1,2), dot(eta_new_i1,xi_2,2), dot(eta_new_i1,nu_vect_inter,2);...
        dot(eta_new_i2,xi_1,2) dot(eta_new_i2,xi_2,2),dot(eta_new_i2,nu_vect_inter,2);...
        dot(sigma_i,xi_1,2),dot(sigma_i,xi_2,2), -dot(sigma_i,nu_vect_inter,2)]; 
        % We used the relation sig_i.nu = - sig_r.nu, sig_i.xi2 = sig_r.xi2
        % and etai2.nu = etar2.nu, etai2.xi2 = -etar2.xi2
        %Reflected beam coordinates with respect to the global cordinates.
        Rot_lbc_eta_sig_r = Rot_etsig_nu_r * Rot_lcs_xi_nu; 
        sigma_r = Rot_lbc_eta_sig_r(3,:); 
        eta_new_r2 = Rot_lbc_eta_sig_r(2,:);
        eta_new_r1 = Rot_lbc_eta_sig_r(1,:);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Reflected GB Plotting in the plane of incidence
         %Plotting the plane of incidence 
%         l_b_crr = 0:l_b_inter/1000:l_b_inter;
%         Eta_i1_path = x_0+transpose(l_b_crr).*eta_i1; %Line equation of all the incident field.
%         plot3(Eta_i1_path(:,1),Eta_i1_path(:,2),Eta_i1_path(:,3),'m','linewidth',5);
%         hold on;
%         l_b_crr = 0:l_b_inter/1000:l_b_inter;
%         X_i_path = x_0+transpose(l_b_crr).*sigma_i; %Line equation of all the incident field.
%         plot3(X_i_path(:,1),X_i_path(:,2),X_i_path(:,3),'r','linewidth',5);
%         hold on;
%         l_b_crr = 0:l_b_inter/1000:l_b_inter;
%         Eta_i2_path = x_0+transpose(l_b_crr).*eta_i2; %Line equation of all the incident field.
%         plot3(Eta_i2_path(:,1),Eta_i2_path(:,2),Eta_i2_path(:,3),'g','linewidth',5);
%         hold on;
        %x = -2000:0.5:2000;
        %y = x; 
        %[X,Y] = meshgrid(x,y);
        %Z = dot(QS_inter,xi_2)/xi_2(3)- (xi_2(1)/xi_2(3))*X- (xi_2(2)/xi_2(3)) *Y;
        %mesh(X,Y,Z);
        %hold on;
%         l_b_crr = 0:l_b_inter/1000:l_b_inter;
%         Eta_i1_new_path = QS_inter +transpose(l_b_crr).*eta_new_i1; %Line equation of all the incident field.
%         plot3(Eta_i1_new_path(:,1),Eta_i1_new_path(:,2),Eta_i1_new_path(:,3),'--m','linewidth',5);
%         hold on;
        %l_b_crr = 0:l_b_inter/1000:l_b_inter;
        %X_i_path = x_0+transpose(l_b_crr).*sigma_i; %Line equation of all the incident field.
        %plot3(X_i_path(:,1),X_i_path(:,2),X_i_path(:,3),'--r','linewidth',2);
        %hold on;
%         l_b_crr = 0:l_b_inter/1000:l_b_inter;
%         Eta_i2_new_path = QS_inter + transpose(l_b_crr).*eta_new_i2; %Line equation of all the incident field.
%         plot3(Eta_i2_new_path(:,1),Eta_i2_new_path(:,2),Eta_i2_new_path(:,3),'--g','linewidth',5);
%         hold on;
%         d = 0:1:l_b_inter;
%         Nu_r_path = QS_inter + transpose(d).*nu_vect_inter;
%         plot3(Nu_r_path(:,1),Nu_r_path(:,2),Nu_r_path(:,3),'--k','linewidth',2);
%         hold on;
%         sig_r_Q = sigma_i - 2*(dot(sigma_i,nu_vect_inter))*nu_vect_inter; %Computing the reflected field direction, the reflected field shall acts as the new incident GB for our algorithm
%         sig_r_Q = sig_r_Q/norm(sig_r_Q);
%         X_r_path = QS_inter+transpose(l_b_crr).*sig_r_Q;
%         plot3(X_r_path(:,1),X_r_path(:,2),X_r_path(:,3),'g','linewidth',2);
%         l_b_crr = 0:l_b_inter/1000:l_b_inter;
%         Xi_1 = QS_inter +transpose(l_b_crr).*xi_1;
%         plot3(Xi_1(:,1),Xi_1(:,2),Xi_1(:,3),'-.b','linewidth',5);
%         hold on;
%         l_b_crr = 0:l_b_inter/1000:l_b_inter;
%         Xi_2 = QS_inter +transpose(l_b_crr).*xi_2;
%         plot3(Xi_2(:,1),Xi_2(:,2),Xi_2(:,3),'-.b','linewidth',5);
%         hold on;
%         l_b_crr = 0:l_b_inter/1000:l_b_inter;
%         Eta_r1_new_path = QS_inter +transpose(l_b_crr).*eta_new_r1; %Line equation of all the incident field.
%         plot3(Eta_r1_new_path(:,1),Eta_r1_new_path(:,2),Eta_r1_new_path(:,3),'-.m','linewidth',5);
%         hold on;
        %l_b_crr = 0:l_b_inter/1000:l_b_inter;
        %X_r_path = QS_inter + transpose(l_b_crr).*sigma_r; %Line equation of all the incident field.
        %plot3(X_r_path(:,1),X_r_path(:,2),X_r_path(:,3),'-.r','linewidth',2);
        %hold on;
%         l_b_crr = 0:l_b_inter/1000:l_b_inter;
%         Eta_r2_new_path = QS_inter +transpose(l_b_crr).*eta_new_r2; %Line equation of all the incident field.
%         plot3(Eta_r2_new_path(:,1),Eta_r2_new_path(:,2),Eta_r2_new_path(:,3),'-.g','linewidth',5);
%         hold on;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % The behaviour of the Gamma(sigma_i) along the sigma_i axis
        % Analysis for the Gamma_r for the reflected beam
        J = [1 0;0 -1]; % Mirror symmetric matrix s
        C = [1/c_radius,0;0,1/c_radius]; % this is only true for a sphere, for a general object it needs to be calculate locally.s
        Gamma_r_temp = J*(Gamma_new_i_Q)*J + 2*cos(theta_i)*J*transpose(inv(Theta_eta_i))*C*(inv(Theta_eta_i))*J;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % New Iteration for the next Interaction for multiple target
        % configuration until the beam emergers out of the target domain
        x_0 = QS_inter;% The Qs will serve as the next x_0 for our new interation in the target domain.
        sig_r_Q = sigma_i - 2*(dot(sigma_i,nu_vect_inter))*nu_vect_inter; %Computing the reflected field direction, the reflected field shall acts as the new incident GB for our algorithm
        sig_r_Q = sig_r_Q/norm(sig_r_Q);
        sigma_i = sig_r_Q; %We shall reiterate our algorithm untill we emerge out of the target domain 
        Gamma_i = Gamma_r_temp; %Relfected Beam's Complex curvature matrix.
        eta_i1 = eta_new_r1; %Transevere coordinate of the reflected field etar1
        eta_i2 = eta_new_r2; %Transevere coordinate of the reflected field etat2
        hit_tot = hit_tot + 1;
    end
end
end