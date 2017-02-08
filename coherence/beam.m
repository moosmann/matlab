ring_energy_petra3 = 6 * 1e9 ; % eV
dist_source_dcm = 50.9;
dist_source_sample = 82.7;
dist_sample_detector = dist_source_sample - dist_source_dcm;

% Lorentz factor
gamma = ring_energy_petra3 / (electron_rest_mass('eV_c2') );

% Magnification factor
mag = (dist_source_sample) / dist_source_dcm;

% Defocusing distance
defocus = dist_sample_detector / mag;

fprintf( '\n%10s:', 'E / GeV')
fprintf( '%10.2g', ring_energy_petra3 / 1e9) 

fprintf( '\n%10s:', 'gamma')
fprintf( '%10.2g', gamma) 

fprintf( '\n%10s:', 'opening angle')
fprintf( '%10.2g', 1/gamma) 

fprintf( '\n%10s:', 'magnification factor')
fprintf( '%10.2g', mag) 

fprintf( '\n%10s:', 'defocusing distance')
fprintf( '%10.2g', defocus) 

fprintf('\n')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lambda = EnergyConverter(30);
k = 2*pi/lambda;

% NIST
Z_A = 0.49954; % atomic-number-to-mass ration
density_nist = 1.7; % g/cm^3 used for muen_rho
mu_rho_30keV = 2.562e-01; % cm^2/g
muen_rho_30keV = 6.614e-02; % cm^2/g

mu_30keV = mu_rho_30keV * density_nist; % 1/cm

fprintf( '\nNIST:\n mu_rho = %g cm^2/g', mu_rho_30keV)
fprintf( '\n mu = %g cm^-1', mu_30keV)

% Henke
density2 = 1.7; % g/cm^3
delta_30keV_d1p7 = 3.91764104e-07;
beta_30keV_d1p7 = 1.21400473e-10;

fprintf( '\nHenke:\n beta = %g ', beta_30keV_d1p7)
fprintf( '\n mu = %g cm^-1', 2 * k * beta_30keV_d1p7 / 100)
%fprintf( '\n beta / density = %g g/cm^3', beta_30keV_d1p7 /density2)



% NIST 2
density = 2.2; % g/cm^-3
en_30keV = 3.001405e+01; % keV
mu_30keV = 5.6916e-01; % cm^-1, mu total
en_1keV = 9.987612e-01;  
mu_1keV = 4.5773E+03; % cm^-1, mu total

% Henke
density1 = 2.2; % g/cm^3
delta_30keV_d2p2 = 5.06988727e-07;
delta_1keV_d2p2 = 0.000481550494;
beta_30keV_d2p2 = 1.57106481e-10;
beta_1keV_d2p2 = 4.8121663e-05;



lambda_30keV = EnergyConverter(30);
lambda_1keV = EnergyConverter(1);
k_30keV = 2*pi/lambda_30keV;
k_1keV = 2*pi/lambda_1keV;

fprintf( '\nNIST : mu(30keV) = %g cm^-1', mu_30keV)
fprintf( '\nHenke: mu(30keV) = %g cm^-1', 2*k_30keV*beta_30keV_d2p2 / 100)
fprintf( '\nNIST : mu(1keV) = %g cm^-1', mu_1keV)
fprintf( '\nHenke: mu(1keV) = %g cm^-1', 2*k_1keV*beta_1keV_d2p2 / 100)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('\n')