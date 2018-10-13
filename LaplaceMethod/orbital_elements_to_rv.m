function [ r_vec, v_vec ] = orbital_elements_to_rv( a, e, i, Omega, w, f  )

mu =  3.986*10^5; %km^3/s^2
p = a*(1-e^2);
h = sqrt(p*mu);

r = p/(1+e*cosd(f));

r_p(1) = r*cosd(f);
r_p(2) = r*sind(f);
r_p(3) = 0;

v_p(1) = -mu/h*sind(f);
v_p(2) = mu/h*(e+cosd(f));
v_p(3) = 0;

DCM_PQW_TO_ECI = [cosd(Omega)*cosd(w)-sind(Omega)*sind(w)*cosd(i), -cosd(Omega)*sind(w)-sind(Omega)*cosd(w)*cosd(i), sind(Omega)*sind(i);
    sind(Omega)*cosd(w)+cosd(Omega)*sind(w)*cosd(i), -sind(Omega)*sind(w)+cosd(Omega)*cosd(w)*cosd(i), -cosd(Omega)*sind(i);
    sind(w)*sind(i), cosd(w)*sind(i), cosd(i)];

r_vec = DCM_PQW_TO_ECI*r_p';
v_vec = DCM_PQW_TO_ECI*v_p';

end