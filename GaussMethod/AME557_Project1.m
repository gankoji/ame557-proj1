
format long g
clear all
clc

%Read observations from file
load OBS.txt
obs.JD   = OBS(:,1);
obs.ra = OBS(:,2);
obs.dec     = OBS(:,3);
clear OBS;

lat = 32.37416; %Observation site latitude (deg);
long = 111.01694; %Observation site longitude (deg);
alt = .757; %Observation site altitude (km)
TOF = 80*1440; %minutes
lst = JD2GMST(obs.JD) -long;
TOF = TOF/60; %translate to hours



output = OrbitCompGauss(lat, lst(1:3), alt, obs.ra(1:3), obs.dec(1:3), obs.JD(1:3), TOF);

output = output';