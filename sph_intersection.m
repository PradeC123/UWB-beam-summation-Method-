% Date - 26/08/2020 (19:00) Israel Daylight Time

% Last update added  "l_b" the distance from the origin of the beam to the
% point of impact on the Sphere
% Date - 30/12/2020 (13:46) Israel Daylight Time (Updates)
%-------------------------------------------------------------------------%
% Sphere  : ||x-x_0||^2 = a^2 ,,, x_0 = (x_0,y_0,z_0) 
% a is the radius of the
% circle 
% Line : x = o + du .
% 1)d is the distance along the line from the point o. 
% 2)u is the direction of the line.
% 3)x is the points on the line
% 4)o is the origin of the line.
% Quadratic equation is ad^2 + bd + c = 0 
% 1) a = u.u   . is the dot product 
% 2) b = 2(u.(o-x_0)) 
% 3) c = (o-x_0).(o-x_0) - a^2 
% D = (u.(o-x_0))^2 - c^2 
% 1) D = 0 (One solution) 2)D > 0 (Two solution) 
% 3) D < 0 (No solution) 
% d = (-2(u.(o-x_0) +-
% sqrt((2(u.(a-x_0)))^2-4u.u(||o-x_0||^2-a^2))/(2||u||^2)

%-------------------------------------------------------------------------%
% INPUT of the function 
% x_0 is the 1X3 matrix, the point where the beam origin.

% sigma_i is the 1X3 matrix, that descrbibes the unit vector in the  direction of the beam propagation

% c is the 1X3 matrix, that describes the center of the sphere.

% r is the 1X1 matrix, that describes the radius of the sphere.
%-------------------------------------------------------------------------%
% OUTPUT of the function.
% ht_pnt is the 1X3 matrix, that describes the point of intersection between a sphere and a line.

% nu_vect is the 1X1 matrix,that describes the unit normal vector in the outward direction to the sphere 

% l_b is the 1X1 matrix, that describes the distance from the origin of the
% beam to the hit point Q on the sphere.

% hit is a 1X1 matrix that descrbies whether the line hits the sphere or
% not
%-------------------------------------------------------------------------%
function [QS,nu_vect,l_b,hit] = sph_intersection(x_0,sigma_i,c,r)
Delta = (dot(sigma_i,x_0-c,2)).^2 - (vecnorm(x_0-c,2,2).^2 -r^2);
if Delta < 0 % Then we find the intersection with the plane at z = z_b + 2*a\
    theta_i = atan2(sigma_i(1),sigma_i(2));
    d = 2*r/(cos(theta_i));
    QS = x_0 + d*sigma_i;
    nu_vect = [0,0,1];
    hit = 0; 
    l_b = d;
else
    %if Delta > 0 %this case corresponds to 2 hits on the sphere, hence we need to consider the first hit.
    d_1 = -dot(sigma_i,x_0-c,2) + sqrt(Delta);
    d_2 = -dot(sigma_i,x_0-c,2) - sqrt(Delta);
    d = min(abs(d_1),abs(d_2)); %the minimum d will correspond to the first hit point on the sphere 
    QS = x_0 + d.*sigma_i;
    nu_vect = (QS - c)./vecnorm((QS - c),2,2);
    l_b = vecnorm(QS-x_0,2,2);    
    hit = not((Delta<0));
    %end
end
end



