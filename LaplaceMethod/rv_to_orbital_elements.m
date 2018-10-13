function [ a, e, i, Omega, w, f ] = rv_to_orbital_elements( r_vec, v_vec )

mu =  3.986*10^5; %km^3/s^2

r = norm(r_vec);
v = norm(v_vec);

Energy = 0.5*v^2 - mu/r;

a = -mu/(2*Energy);  %km

dot_r_v = dot(r_vec,v_vec);

e_vec = 1/mu*((v^2-mu/r)*r_vec - dot_r_v*v_vec);

e = norm(e_vec);

h_vec = cross(r_vec,v_vec);

h = norm(h_vec);

i = acosd(h_vec(3)/h);    %deg

n_vec = cross([0 0 1],h_vec);

n = norm(n_vec);

if e == 0
    Omega = [];
else
    cos_Omega = n_vec(1)/n;
    sin_Omega = n_vec(2)/n;
    Omega = atan2d(sin_Omega,cos_Omega);   %deg
    
    if Omega < 0
        Omega = Omega + 360;
    end
end

if e_vec(3) > 0
    w = acosd(dot(n_vec,e_vec)/(n*e));   %deg
elseif e_vec(3) < 0
    w = 360 - acosd(dot(n_vec,e_vec)/(n*e));   %deg
end

if i == 0 && e > 0
    w = Omega + w;   %deg
elseif e == 0
    w = [];   %deg
end

if dot_r_v > 0
    f = acosd(dot(e_vec,r_vec)/(e*r));   %deg
elseif dot_r_v < 0
    f = 360 - acosd(dot(e_vec,r_vec)/(e*r));   %deg
end

if i > 0 && e == 0
    f =  w + f;   %deg
elseif i == 0 && e == 0
    f = Omega + w + f;   %deg
end
end