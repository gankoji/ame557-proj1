function output = Kepler_Prob(a,e,f1,dt_hours)


dt= dt_hours*60*60; %sec
mu = 3.986*10^5;

f1_rad = degtorad(f1);

n = sqrt(mu/a^3)
E1 = 2*atan(sqrt((1-e)/(1+e))*tan(f1_rad/2))
M1 = E1-e*sin(E1)
M2 = M1+n*dt

E2_0 = M2-e*sin(M2)+e^2/2*sin(2*M2)

err = E2_0 -e*sin(E2_0)-M2
while abs(err)>10^-9
    disp('i')
    df = 1- e*cos(E2_0);
    E2 = E2_0-err/df;
    err = E2-e*sin(E2) - M2
    E2_0 = E2;
    
end
f2 = 2*atan(sqrt((1+e)/(1-e))*tan(E2/2))

r2 =  a*(1-e^2)/(1+e*cos(f2))
gamma = atan((e*sin(f2))/(1+e*cos(f2)));
E = mu/(2*a)
v2 = sqrt(2*(E+mu/r2));

output = [radtodeg(f2)];

