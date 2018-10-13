function output = RVtoOE(r_vec,v_vec)

mu = 3.986*10^5;
r = norm(r_vec);
v = norm(v_vec);

%Energy
E = v^2/2 - mu/r;

%Semimajor axis
a = -mu/(2*E)

%Eccentricity
e_vec = 1/mu*((v^2-mu/r)*r_vec-(r_vec'*v_vec)*v_vec);
e = norm(e_vec);

%Momentum
h = cross(r_vec,v_vec);
%Inclination
i = acosd(h(3)/norm(h));

%Longitude of the ascending node
n = cross([0;0;1],h);
if i ==0
    Omega = [];
else
    Omega = acosd(([1,0,0]*n)/norm(n));

    if n(2) < 0 
        Omega = 360 - Omega;
    end
end

%Arguement of periapsis
if norm(n)==0 && e ~= 0; %elliptical equatorial orbit
    I = [1;0;0];
    w = acosd((I'*e_vec)/(norm(I)*e));
else 
    w = acosd((n'*e_vec)/(norm(n)*e));

    if e_vec(3) < 0 
        w = 360 - w;
    end
end

%True anomaly at epoch
if e ~= 0
    nu = acosd((e_vec'*r_vec)/(e*r));
    if (dot(r_vec,v_vec))<0
        nu = 360 - nu;
    end

end

output = [a,e,i,Omega,w,nu];
end