% formulas: https://physics.nist.gov/PhysRefData/XrayMassCoef/chap2.html
% values: https://physics.nist.gov/PhysRefData/XrayMassCoef/ComTab/bone.html
% constants: https://physics.nist.gov/PhysRefData/XrayMassCoef/tab2.html

% Energy / MeV
E = [1.00000E-02 1.50000E-02 2.00000E-02 3.00000E-02 4.00000E-02 5.00000E-02 6.00000E-02 8.00000E-02];

% X-ray mass attenuation coefficient mu/rho / (cm^2/g)
mu_over_rho = [ 2.851E+01 9.032E+00 4.001E+00 1.331E+00 6.655E-01 4.242E-01 3.148E-01 2.229E-01];

% transmission I/I_0 = exp( - mu * t )
t = 0.2; % thickness / cm
rho = 1.920E+00; % density / (g/cm^3)
trans = exp( - mu_over_rho * rho * t );

%% Figure
x = 1:5;
figure( 'Name', 'transmission of cortical bone' )
plot( 1000 * E(x), trans(x) )
xlabel( 'keV' )
ylabel( 'tranmission' )

