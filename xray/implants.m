% Linear extrapolation of transmission of Mg10Gd screw at 40 keV
% using transmission of a solid at 1 bar ~ 750 Torr and Mg10Gd where 
% 10 refers to 10 mass-percent. With a standard atomic weight of 
u_Mg = 24.3;
% and 
u_Gd = 157.3; 
%, the relative amount of substance is:
% m_Gd = 0.1 m_Mg
% N_Gd * u_Gd = 0.1 * N_Mg * u_Mg
% With N_Gd = 1:
N_Mg = 10 * u_Gd / u_Mg;
disp(N_Mg)
% the stochiometric number is Mg64.7Gd1

% Transmission values at two nearby points below 30 keV from Henke
e1=25000;
t1=0.75525E-01;

e2=30000;
t2=0.20745;

% Linear slope
m=(t2-t1)/(e2-e1);

% Transmission at 40 keV
t3 = t1 + m*(40000-e1);
