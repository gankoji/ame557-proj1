function [ECI_vec] = Geodetic2ECI(lat, alt, lst)

R_E     = 6378.145; % km (Equatorial radius)
e       = 0.08182;
f       = 0.0033528;

N = R_E/sqrt(1-(2*f-f^2)*sind(lat)^2);  % radius of curvature in the prime vertical

X = (N+alt)*cosd(lat)*cosd(lst);
Y = (N+alt)*cosd(lat)*sind(lst);
Z = (N*(1-f)^2+alt)*sind(lat);

ECI_vec = [X;Y;Z];
end