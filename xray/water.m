% formulas:
% values:
% constants: https://physics.nist.gov/PhysRefData/XrayMassCoef/tab2.html

% Energy / MeV
E = [1.00000E-03 1.50000E-03 2.00000E-03 3.00000E-03 4.00000E-03 5.00000E-03 6.00000E-03 8.00000E-03 1.00000E-02 1.50000E-02 2.00000E-02 3.00000E-02 4.00000E-02 5.00000E-02 6.00000E-02 8.00000E-02 1.00000E-01 1.50000E-01 2.00000E-01 3.00000E-01 4.00000E-01 5.00000E-01 6.00000E-01 8.00000E-01 1.00000E+00 1.25000E+00 1.50000E+00 2.00000E+00 3.00000E+00 4.00000E+00 5.00000E+00 6.00000E+00 8.00000E+00 1.00000E+01 1.50000E+01 2.00000E+01];
% X-ray mass attenuation coefficient mu/rho / (cm^2/g)
mu_over_rho = [4.078E+03 1.376E+03 6.173E+02 1.929E+02 8.278E+01 4.258E+01 2.464E+01 1.037E+01 5.329E+00 1.673E+00 8.096E-01 3.756E-01 2.683E-01 2.269E-01 2.059E-01 1.837E-01 1.707E-01 1.505E-01 1.370E-01 1.186E-01 1.061E-01 9.687E-02 8.956E-02 7.865E-02 7.072E-02 6.323E-02 5.754E-02 4.942E-02 3.969E-02 3.403E-02 3.031E-02 2.770E-02 2.429E-02 2.219E-02 1.941E-02 1.813E-02];

% transmission I/I_0 = exp( - mu * t )
t = 6.5; % thickness / cm
rho = 1E+00; % density / (g/cm^3)
trans = exp( - mu_over_rho * rho * t );

%% Figure
x = 1:numel(E);
x = 13:15;
figure( 'Name', 'transmission of water' )
plot( 1000 * E(x), trans(x) )
xlabel( 'keV' )
ylabel( 'tranmission' )

%% discrete values similiar to Henke
x = [7 8 9 11 12];
E_eV = E(x) * 1e6;
disp(  E_eV );
disp( trans(x) )