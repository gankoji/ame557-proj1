function [ r2, v2, gamma2, f2 ] = propagate_position( a, e, f, flight_time )

% AME557, HW2, Problem 4.16
% Gustavo Lee
mu = 3.986*10^5;    %km^3/s^2
n = sqrt(mu/a^3);

E1 = 2*atan2(sqrt(1-e)*tand(f/2),sqrt(1+e));

M1 = E1 - e*sin(E1);

M2 = M1 + n*flight_time;

% Newton's method
E(1) = M2;
i = 1; err = 9999;
while abs(err) > 0.001
    E(i+1) = E(i) - (E(i)-e*sin(E(i))-M2)/(1-e*cos(E(i)));
    err = E(i+1)-E(i);
    i = i+1;
end
E2 = E(i);

f2 = 2*atan2d(sqrt(1+e)*tan(E2/2),sqrt(1-e));      %deg

r2 = a*(1-e^2)/(1+e*cosd(f2));      %km

Energy = -mu/(2*a);

v2 = sqrt(2*(Energy+mu/r2));      %km/s

gamma2 = atan2d(e*sind(f2),1+e*cosd(f2));    %deg

end