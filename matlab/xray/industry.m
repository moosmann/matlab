
eV_to_J = 1 / 6.24e18;

E = (50:1:200)*1e3;
  
% Energy in eV
% X-ray mass attenuation coefficient mu/rho in m^2 / kg
[energy_hf, ~, mac_hf] = read_nist_txt( 'hafnium' );
[energy_ag, ~, mac_ag] = read_nist_txt( 'silver' );


% thickness in m
thickness_hf = ((2)*1e-3)'; 
thickness_ag = ((2)*1e-3)'; 

% Hf density in kg / m^3
density_hf = 1.331E+01 * 1000;
density_ag = 1.050E+01 * 1000;

% transmission I/I_0 = exp( - mu * t )
transmission_hf = exp( - interp1( energy_hf, mac_hf, E ) * density_hf .* thickness_hf );
transmission_ag = exp( - interp1( energy_ag, mac_ag, E ) * density_ag .* thickness_ag );

%fprintf( '\ntransmission_water : %f', transmission_hf )

% subplot(1,3,1)
% plot(E / 1e3, transmission_hf)
% xlabel('E / keV')
% ylabel('transmission')
% title('Hafnium')
% 
% subplot(1,3,2)
% plot(E / 1e3, transmission_ag)
% xlabel('E / keV')
% ylabel('transmission')
% title('Silver')
% 
% subplot(1,3,3)
% plot(E / 1e3, transmission_ag .* transmission_hf)
% xlabel('E / keV')
% ylabel('transmission')
% title('Silver * Hafnium')
% 

Y = [transmission_ag; transmission_hf; transmission_ag .* transmission_hf];
plot(E / 1e3, Y )
xlabel('E / keV')
ylabel('transmission')
%title('Hafnium')
legend( {'Silver', 'Hafnium', 'Hf * AG'} )
grid on
