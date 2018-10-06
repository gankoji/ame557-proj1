%[r0, v0, oe0, rf, vf, oef] = 
function output = OrbitCompGauss(lat, lst, alt, ra, dec, JD, TOF)
%Inputs
%lat - Latitude of the observation site (in degrees)
%lst - Local siderial time of observation site (in degrees)
%alt - Altitude above sea level of the observation site
%ra - 3x1 vector of Right ascensions in SEZ frame (in degrees)
%dec - 3x1 vector of declinations in SEZ frame (in degrees)
%JD - 3x1 vector of Julian Dates corresponding to the 3 observations
%TOF - Time of flight for orbit propagation

mu =  3.986*10^5; %km^3/s^2
Re = 6378.137; %km
day2sec = 24*60*60;
tau1 = (JD(1) - JD(2))*day2sec %seconds
tau3 = (JD(3) - JD(2))*day2sec %seconds

a1 = tau3/(tau3-tau1)
a3 = -tau1/(tau3-tau1)
a1u = (tau3*((tau3-tau1)^2 - tau3^2))/(6*(tau3-tau1)) %seconds^2
a3u = -(tau1*((tau3-tau1)^2 - tau1^2))/(6*tau3-tau1) %seconds^2

L = [cosd(dec').*cosd(ra'); cosd(dec').*sind(ra');sind(dec')]

r_site = (Re+alt)*[cosd(lat)*cosd(lst');cosd(lat)*sind(lst');sind([lat,lat,lat])]

M = L\r_site


d1 = M(2,1)*a1-M(2,2)+M(2,3)*a3 %km
d2 = M(2,1)*a1u + M(2,3)*a3u %km*sec^2

C = L(:,2)'*r_site(:,2)
 
%Coefficients of eigth order polynomial
r2_poly = [1 0 -(d1^2+2*C*d1+norm(r_site(:,2))^2) 0 0 -2*mu*(C*d2+d1*d2) 0 0 -mu^2*d2^2] 

sol = roots(r2_poly)
r2_sol = sol(imag(sol)==0 & sol>=0) %selects the positive real solution

u = mu/r2_sol^3 %sec^-2

c_vec = -[-a1-a1u*u;1;-a3-a3u*u]

p = -(M*c_vec)./c_vec

r1 = p(1)*L(:,1)+r_site(:,1)
r2 = p(2)*L(:,2)+r_site(:,2)
r3 = p(3)*L(:,3)+r_site(:,3)

t31 = JD(3)-JD(1);
t32 = JD(3)-JD(2);
t21 = JD(2)-JD(1);

Z23 = cross(r2,r3)

a_cop = 90 - acosd((Z23'*r1)/(norm(Z23)*norm(r1)));

v2 = -t32*(1/(t21*t31)+mu/(12*norm(r1)^3))*r1 + (t32-t21)*(1/(t21*t32)
    +mu/(12*norm(r2)^3))*r2 + t21*(1/(t32*t31)+mu/(12*norm(r3)^3))*r3;


output = 0;





