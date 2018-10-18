% formulas: https://physics.nist.gov/PhysRefData/XrayMassCoef/chap2.html
% values: https://physics.nist.gov/PhysRefData/XrayMassCoef/ComTab/air.html
% constants: https://physics.nist.gov/PhysRefData/XrayMassCoef/tab2.html

% Energy / MeV
E = [
1.00000E-03  
1.50000E-03  
2.00000E-03  
3.00000E-03  
3.20290E-03  
3.20290E-03  
4.00000E-03  
5.00000E-03  
6.00000E-03  
8.00000E-03  
1.00000E-02  
1.50000E-02  
2.00000E-02  
3.00000E-02  
4.00000E-02  
5.00000E-02  
6.00000E-02  
8.00000E-02  
1.00000E-01  
1.50000E-01  
2.00000E-01  
3.00000E-01  
4.00000E-01  
5.00000E-01  
6.00000E-01  
8.00000E-01  
1.00000E+00  
1.25000E+00  
1.50000E+00  
2.00000E+00  
3.00000E+00  
4.00000E+00  
5.00000E+00  
6.00000E+00  
8.00000E+00  
1.00000E+01  
1.50000E+01  
2.00000E+01  
    ];
E = E(:)';

% X-ray mass attenuation coefficient mu/rho / (cm^2/g)
mu_over_rho = [
3.606E+03   
1.191E+03   
5.279E+02   
1.625E+02   
1.340E+02   
1.485E+02   
7.788E+01   
4.027E+01   
2.341E+01   
9.921E+00   
5.120E+00   
1.614E+00   
7.779E-01   
3.538E-01   
2.485E-01   
2.080E-01   
1.875E-01   
1.662E-01   
1.541E-01   
1.356E-01   
1.233E-01   
1.067E-01   
9.549E-02   
8.712E-02   
8.055E-02   
7.074E-02   
6.358E-02   
5.687E-02   
5.175E-02   
4.447E-02   
3.581E-02   
3.079E-02   
2.751E-02   
2.522E-02   
2.225E-02   
2.045E-02   
1.810E-02   
1.705E-02   
    ];
mu_over_rho = mu_over_rho(:)';

% transmission I/I_0 = exp( - mu * t )
t = 1; % thickness / cm
rho = 1.205E-03; % density / (g/cm^3)
trans = exp( - mu_over_rho * rho * t );

%% Figure
figure( 'Name', 'transmission of air' )
plot( 1000 * E, trans )
xlabel( 'keV' )
ylabel( 'tranmission' )

%% discrete values similiar to Henke
x = 1:14;
E2 = E(x) * 1e3;
trans2 = trans(x);
disp(  E2 );
disp( trans2 )
figure( 'Name', 'transmission of air in the range of Henke' )
plot( E2, trans2 )