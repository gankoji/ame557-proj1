function [r0,v0,oe0, rf, vf, oef] = OrbitCompLaplace(lat, lst, alt, ra, dec, JD, TOF)
%Inputs
%lat - Latitude of the observation site (in degrees)
%lst - Local sidereal time of observation site (in degrees)
%alt - Altitude above sea level of the observation site
%ra - 3x1 vector of Right ascensions in SEZ frame (in degrees)
%dec - 3x1 vector of declinations in SEZ frame (in degrees)
%JD - 3x1 vector of Julian Dates corresponding to the 3 observations
%TOF - Time of flight for orbit propagation (minutes)

for i = 1:numel(JD)
    L(:,i) = [cosd(dec(i))*cosd(ra(i));...
        cosd(dec(i))*sind(ra(i));...
        sind(dec(i))];
end
for i = 1:numel(JD)
    Ldot(:,i) = ((2*JD(i)-JD(2)-JD(3))/((JD(1)-JD(2))*(JD(1)-JD(3)))*L(:,1) + ...
        (2*JD(i)-JD(1)-JD(3))/((JD(2)-JD(1))*(JD(2)-JD(3)))*L(:,2) + ...
        (2*JD(i)-JD(1)-JD(2))/((JD(3)-JD(1))*(JD(3)-JD(2)))*L(:,3))/86400;
    
    Ldotdot(:,i) =  (2/((JD(1)-JD(2))*(JD(1)-JD(3)))*L(:,1) + ...
        2/((JD(2)-JD(1))*(JD(2)-JD(3)))*L(:,2) + ...
        2/((JD(3)-JD(1))*(JD(3)-JD(2)))*L(:,3))/(86400*86400);
end

mu =  3.986*10^5; %km^3/s^2
OMEGA_E  = 7.2921150E-5; % rad/s
OMEGA_Evec = [0;0;OMEGA_E];

ECI_vec       = Geodetic2ECI(lat, alt, lst(2));
ECI_vecdot    = cross(OMEGA_Evec,ECI_vec);
ECI_vecdotdot = cross(OMEGA_Evec,ECI_vecdot);

R = norm(ECI_vec);
L_dot_R = dot(L(:,2), ECI_vec);
D  = 2*det([L(:,2),Ldot(:,2), Ldotdot(:,2)]);
D1 = det([L(:,2),Ldot(:,2),ECI_vecdotdot]);
D2 = det([L(:,2),Ldot(:,2),ECI_vec]);
D3 = det([L(:,2),ECI_vecdotdot,Ldotdot(:,2)]);
D4 = det([L(:,2),ECI_vec,Ldotdot(:,2)]);

c_0 = -4*mu^2*D2^2;
c_3 = -8*mu*D1*D2+4*mu*D*D2*L_dot_R;
c_6 = -4*D1^2+4*D*D1*L_dot_R-R^2*D^2;
c_8 = D^2;

Z = roots([c_8 0 c_6 0 0 c_3 0 0 c_0]);

for i=1:numel(Z)
    if (isreal(Z(i)) && Z(i) > 0)
        r = Z(i);
        break;
    end
end

rho = -2*D1/D-2*mu*D2/(r^3*D);
rhodot = -D3/D - mu/r^3*(D4/D);

r0 = rho*L(:,2)+ECI_vec;
v0 = rhodot*L(:,2)+rho*Ldot(:,2)+ECI_vecdot;


[ a, e, i, Omega, w, f ] = rv_to_orbital_elements(r0,v0);
oe0=[a;e;i;Omega;w;f];

% Convert TOF to seconds for propagation
[ ~, ~, ~, f2 ] = propagate_position( a, e, f, TOF*60 );

[ rf, vf ] = orbital_elements_to_rv( a, e, i, Omega, w, f2  );
oef=[a;e;i;Omega;w;f2];