% sample consists mainly of water as well as the cured agaros solution, container PMMA or similiar stuff
% container thickness = dx ~ 12 mm = 12000 microns
% container wall thickness ~ 0.5 mm
dx = 0.012; % m

% Energy:
% h      = 6.62606896e-34;
% c      = 299792458;
% q      = 1.60217733e-16;
% lambda = h*c./(E/keV*q)
% k = 2*pi/lambda
E22 = 22;
lambda22 = EnergyConverter(E22);
k22 = 2*pi/lambda22;
E30 = 30;
lambda30 = EnergyConverter(30);
k30 = 2*pi/lambda30;

% phase shift approximately:
% phi = -k*delta_k*dx
% absorption approximately (Beer's law of absorption)
% mu_k = 2*k*beta_k 
% Filter Transmission = abs(T)
% T = exp(-mu*dx)*exp(phi)
% I_z=0 = abs(T)*I0
% I0 = incident beam = flat field

% from http://henke.lbl.gov:
% H2O Density=1. Thickness=12000. microns
%  Photon Energy (eV), Transmission
%     22000.      0.50867    
%     30000.      0.67819 
% T = exp(-t)
T22 = 0.50867;
t22 = -log(T22);
T30 = 0.67819;
t30 = -log(T30);
% from http://henke.lbl.gov
%  H2O Density=1.
%  Energy(eV), Delta, Beta
%   22000.  4.76455682E-07  2.52630611E-10
%   30000.  2.56114134E-07  1.06431752E-10
delta22 = 4.76455682E-07;
beta22 = 2.52630611E-10;
delta30 = 2.56114134E-07;  
beta30 = 1.06431752E-10;
% Values of T22 is backed up by experimental data. Checking a
% flat-and-dark-field corrected and hot-pixel-filtered projection. A
% central region in the projection where the flats are not affected by the
% flat field is considered and the mean is computed over the region: ~0.52
% for E = 22 keV

% Going from 22 keV to 30 keV: E22 -> E33
% Estimate decrease in counts:
% Using scintillator GGG the absorption effiency drops to 67%
dAE = 0.67;
% The absorption in the samples drops to 75%
dT = T30/T22;
% From energy spectrum of BM05: decrease in photon intensity
dP = 133.8/151.9;
dCounts = dAE*dT*dP;

% Estimate decrease in intensity:
% Assume validity of linearized TIE: I =
% I_z-0*(1-lambda*z*laplace(dphi)), thus intensity contrast g = I_z/I_z=0-1 ~
% lambda*z*phi.
% phi = -k*delta_k*dx
dIC = lambda30/lambda22*k30/k22*delta30/delta22;

% From data:
z1 = 0.123;
z2 = 0.423;
dz = z2/z1;
ICmaxVar1 = 0.116704;
ICmaxVar2 = 0.265036;
dICmaxVar = ICmaxVar2/ICmaxVar1;
ICstd1 = 0.0114892;
ICstd2 = 0.0286326;
dICstd = ICstd2/ICstd1;
% E = 22, pixel size = 1.5e-6. Read flat-field- and hot-pixel-corrected
% projection. Crop central part.
%Domain of int100c: [0.456902,0.573606], Mean=0.515346, Sta=0.0114892, Max-Min=0.116704
%Domain of int400c: [0.379215,0.644251], Mean=0.508195, Sta=0.0286326, Max-Min=0.265036

