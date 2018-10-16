format long g
clear all
clc

% %Read observations from file
% load OBS.txt
% obs.JD  = OBS(:,1);
% obs.ra  = OBS(:,2);
% obs.dec = OBS(:,3);
% clear OBS;

%Read observations from file
load Generated_OBS.txt
obs.JD  = Generated_OBS(:,1)
obs.ra  = Generated_OBS(:,2)
obs.dec = Generated_OBS(:,3)
clear Generated_OBS;

%Site location
lat = 32 + 22/60 + 27/3600;       %Observation site latitude (deg);
long = 360-(111 + 1/60 + 1/3600); %Observation site east longitude (deg);
alt = 757/1000;                   %Observation site altitude (km)
TOF = 80*24*60;                   %TOF (minutes)

lst = JD2GMST(obs.JD) + long;      %Local Sidereal Time


% output = OrbitCompGauss(lat, lst(1:3), alt, obs.ra(1:3), obs.dec(1:3), obs.JD(1:3), 0);
[r0,v0,oe0, rf, vf, oef] = OrbitCompLaplace(lat, lst(1:2:5), alt, obs.ra(1:2:5), obs.dec(1:2:5), obs.JD(1:2:5), TOF)