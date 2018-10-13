function output = OEtoRV(oe0,f)

mu = 3.986*10^5;
a = oe0(1);
e = oe0(2);
i = oe0(3);
omega = oe0(4);
w = oe0(5);

p = a*(1-e^2)
h = sqrt(p*mu)
r0 = p/(1+e*cosd(f));

rpqw = [r0*cosd(f);r0*sind(f);0];

vpqw = [-mu/h*sind(f);mu/h*(e+cosd(f));0];

R1 = [cosd(w), sind(w),0;-sind(w), cosd(w),0;0,0,1];
R2 = [1,0,0;0,cosd(i),sind(i);0,-sind(i),cosd(i)];
R3 = [cosd(omega), sind(omega),0;-sind(omega), cosd(omega),0;0,0,1];

R = (R1*R2*R3)';

r_eci = R*rpqw;
v_eci = R*vpqw;

output = [r_eci;v_eci];

end