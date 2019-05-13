%function out = dose_synchroload(fluxDensity_ph_per_s_mm2,energy.value,absorption_length.value,scan_time.value,object_length.value,objectDensity_g_per_ccm,container_length.value)
% One gray is the absorption of one joule of energy, in the form of
% ionizing radiation, per kilogram of matter.
%
%   1 Gy = 1 J / kg = m^2 / s^2
clear all

%% ID 19
id19.source = 'undulator U13';
id19.spectrum = 'pink';
id19.bandwidth = 0.04;
id19.ring_current.value = 193;
id19.ring_current.unit = 'mA';
id19.peak_energy.value = 26.3;
id19.peak_energy.unit = 'keV';
id19.flux.value = 166e-12;
id19.flux.unit = 'photons / s / mm^2';

%% 2-BM-B
% @ 30keV DMM: Henke: 8.0e13 ph/s/mrad^2/0.1%BW
% flux density = 8e13*10 ph/s/1%BW / (35m)^2 /mrad^2
% = 6.54e11 ph/s/1%BW/mm^2 ~ 10^12 ph/s/mm^2
fluxDensity_ph_per_s_mm2 = 10^12; % photons / s / mm ^2

energy.value = 30 ;
energy.unit = 'keV' ;

% for water: l_att = 1 / (4*pi*beta/lambda) = 3.10 cm at 30 keV
absorption_length.value = 30.9;
absorption_length.unit = 'cm';

scan_time.value = 20 ;
scan_time.unit = 's' ;

object_length.value = 1.0 ;
object_length.unit = 'mm';

object_density.value  = 0.997 ; 
object_density.unit  = 'g / ccm';
object_density.unit_equivalent = '1000 kg / m^3' ;

container_length.value = 12 ; % mm
container_length.unit = 'mm';

%% P05 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% absorber

% CVD 300 micron
% Cu 50 micron
% transmission at 30000 eV: 0.59

% glassy / vitreous carbon / Glaskohlenstoff

absorber(1).type = 'GC';
% density
absorber(1).density.value = [1.4 1.5];
absorber(1).density.unit = 'g / ccm';
% thickness
absorber(1).thickness.value = 0.004;
absorber(1).thickness.unit = 'm';

% energy
p05.energy.value = [14400 30000 ];
p05.energy.unit = 'eV';

% max flux at sample position
p05.flux.value = [7e13, 525542981736 / 3 ]; 
p05.flux.unit = 'photons / s';

% spot size at sample position
p05.area.value = [6.000 * 3.0, 1.0 * 1.0];
p05.area.unit = 'mm^2';

% flux density
p05.flux_density.value = p05.flux.value ./ p05.area.value;
p05.flux_density.unit = 'photons / s / mm^2';

%% Materials

% sample environment

transmission.PEEK = 0.85; %extrapolated guess from 30

% average bone+PEEK
transmission.data.bone_PEEK = 0.61;

%% Bone

% formulas: https://physics.nist.gov/PhysRefData/XrayMassCoef/chap2.html
% constants: https://physics.nist.gov/PhysRefData/XrayMassCoef/tab2.html

bone.thickness.value = 0.001;
bone.thickness.unit = 'm';
bone.density.value = 1.920E+00 / 1000 / 0.01^3 ;
bone.density.unit = 'kg / m^3';
bone.density.unit_equiv = 'g / cm^3 * 1000 * 0.01^3';

% Energy / MeV
energy.value = [1.00000E-02 1.50000E-02 2.00000E-02 3.00000E-02 4.00000E-02 5.00000E-02 6.00000E-02 8.00000E-02];
energy.unit = 'MeV';

% X-ray mass attenuation coefficient mu/rho / (cm^2/g)
mass_att_coeff.value = [ 2.851E+01 9.032E+00 4.001E+00 1.331E+00 6.655E-01 4.242E-01 3.148E-01 2.229E-01];
mass_att_coeff.unit = 'cm^2 / g';
mass_att_coeff.def = 'mu / rho';

% transmission I/I_0 = exp( - mu * t )
trans = exp( - mass_att_coeff.value * bone.density.value * bone.thickness.value * 100 );
mu_over_rho_35000eV = mean( mass_att_coeff.value(4:5) );
transmission_bone_35000eV = exp( - mu_over_rho_35000eV * bone.density.value * bone.thickness.value*100 );
bone.volume.value = (bone.thickness.value)^3;
bone.volume.unit = 'm^3';
bone.mass.value = bone.density.value * bone.volume.value; 
bone.mass.unit = 'kg';

absorption = 1 - transmission_bone_35000eV;

% exposure
exposure_time.value = [1200 * 0.15, 3000 * 0.1];
exposure_time.unit =  's';

energy.value = 34000;
energy.unit = 'eV';

%% Dose
flux_density_30000eV = 1;
dose = transmission.PEEK * absorption * flux_density_30000eV * (bone.thickness.value*1000)^2 * exposure_time.value * energy * eV_to_J / bone.mass.value;

fprintf( '\ndose : %g\b kGy', dose / 1000 )

% Attenuation
out.containerTransmission = exp( - (container_length.value -object_length.value)/2 / absorption_length.value) ;
out.objectAbsorption = 1-exp( - object_length.value/ absorption_length.value) ;
% Number of photons
out.absorbedPhotons = out.containerTransmission * out.objectAbsorption * fluxDensity_ph_per_s_mm2 * scan_time.value  * (object_length.value)^2;
% Cubic mass of object
out.cubeMass_g = object_length.value^3 * objectDensity_g_per_ccm / 1000 ;
% Dose
out.dose_Gy = out.absorbedPhotons * energy.value*10^3 * eV_to_J ...
    / ( out.cubeMass_g / 1000 );

