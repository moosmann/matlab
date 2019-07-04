% One gray is the absorption of one joule of energy, in the form of
% ionizing radiation, per kilogram of matter.
%
% 1 Gy = 1 J / kg = m^2 / s^2

eV_to_J = 1 / 6.24e18;

bin = 2;
energy = 20000:500:50000;
fprintf( ' energy: min : %g keV, max : %g keV, steps : %u', min( energy )/1000, max( energy )/1000, numel( energy ) )
voxel_size.value = bin * 2 * 3.115e-6;
voxel_size.unit = 'm';
fprintf( '\n voxel_size : %g %s\n', voxel_size.value, voxel_size.unit )

exp_time_per_image.value = 300e-3;
exp_time_per_image.unit = 's';

%% Scintillator

cdwo100 = cadmium_tungstate(energy,100e-6);
cdwo300 = cadmium_tungstate(energy,300e-6);
Y = [cdwo100.absorption.value; cdwo300.absorption.value];
if exist( 'h1' , 'var' ) && isvalid( h1 )
    figure(h1)
else
    h1 = figure( 'Name', 'Absorption Scintillator' );
end
plot( energy / 1000, Y)
legend( {'100 micron', '300 micron'} )
xlabel( 'energy / keV' )

%% Bone sample

scan_path = '/asap3/petra3/gpfs/p05/2017/data/11003440/processed/syn32_99R_Mg10Gd_4w/segmentation/2017_11003440_syn32_99R_Mg10Gd_4w_labels';
if ~exist( 'vol', 'var' )
    vol = read_images_to_stack( scan_path );
    domain( vol(:), 1, 'original volume  ' )
    vol = single( vol );
    vol = Binning( vol, bin ) / bin^3;
    domain( vol(:), 1, 'single conversion' )
    vol = permute( vol, [1 2 3] );
end
fprintf( '\nvolume shape : %u %u %u', size( vol ) )
vol_bone = vol == 205;
sc =  (vol == 105 | vol == 5);

bc = bone_cortical(energy,2e-3);

% volume
bc.total_volume.value = sum( vol_bone(:) ) * voxel_size.value^3;
bc.total_volume.unit = 'm^3';
bc.total_mass.value = bc.density.value * bc.total_volume.value;
bc.total_mass.unit = 'kg';

total_volume = voxel_size.value^3 * numel( vol_bone );
bone_volume = bc.total_volume.value;
bone_mass = bc.total_mass.value;
fprintf( '\n bone volume : %g m^3 (%.2g%%)', bone_volume, 100 * bone_volume / total_volume )
fprintf( '\n bone mass : %g mg', bone_mass * 1e6 )

%% Dose
% Projected bone thickness
projected_thickness_b = permute( astra_make_sino_3D( vol_bone ), [ 3 1 2] ) * voxel_size.value;
fprintf( '\n projections shape : %u %u %u', size( absorbed_energy ) )
num_proj = size( projected_thickness_b, 3 );

% Scan time
total_scan_time = num_proj * exp_time_per_image.value;
fprintf( '\n total scan time : %g s', total_scan_time )

% Preallocation
absorbed_energy_per_image = zeros( [1, num_proj] );
absorbed_energy_per_scan = zeros( [1, numel( energy )] );

% Flux values from previous experiment
flux_density.value = 2.5655e+17;
flux_density.unit = 'photons / s / m^2';
area = size( projected_thickness_b, 1) * size( projected_thickness_b, 2) * voxel_size.value^2;
flux.value = flux_density.value * area;
flux.unit = 'photons / s';

fprintf( '\n Start loop over energy and projections' )
bone_density = bc.density.value;
mac = bc.mass_att_coeff.value;
exp_time = exp_time_per_image.value;
flu = flux.value;
fprintf( '\n energy step (%u): ', numel( energy ) )
for kk = 1:numel( energy )
    fprintf( ' %u', kk)
    mac_kk = mac(kk);
    energy_kk = energy(kk);
    for ll = 1:num_proj
        % absorption image
        a = 1 - exp( - bone_density * projected_thickness_b(:,:,ll) * mac_kk );
        % absorbed energy of full image
        absorbed_energy_per_image(ll) = sum( sum( a * flu * exp_time .* energy_kk * eV_to_J ) );
    end
    % absorbed energy of full tomo scan summing over all images
    absorbed_energy_per_scan(kk) = sum( absorbed_energy_per_image );
end
fprintf( '\n Loop finished' )

dose.value = absorbed_energy_per_scan / bone_mass;
dose.unit = 'Gy';

if exist( 'h2' , 'var' ) && isvalid( h2 )
    figure(h2)
else
    h2 = figure( 'Name', 'Dose: Cortical bone' );
end
plot( energy / 1000, dose.value / 1000 )
xlabel( 'energy / keV' )

fprintf( '\n' )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% ID 19
% id19.source = 'undulator U13';
% id19.spectrum = 'pink';
% id19.bandwidth = 0.04;
% id19.ring_current.value = 193;
% id19.ring_current.unit = 'mA';
% id19.peak_energy.value = 26.3;
% id19.peak_energy.unit = 'keV';
% id19.flux.value = 166e-12;
% id19.flux.unit = 'photons / s / mm^2';
%
% %% 2-BM-B
% % @ 30keV DMM: Henke: 8.0e13 ph/s/mrad^2/0.1%BW
% % flux density = 8e13*10 ph/s/1%BW / (35m)^2 /mrad^2
% % = 6.54e11 ph/s/1%BW/mm^2 ~ 10^12 ph/s/mm^2
% fluxDensity_ph_per_s_mm2 = 10^12; % photons / s / mm ^2
%
% energy.value = 30 ;
% energy.unit = 'keV' ;
%
% % for water: l_att = 1 / (4*pi*beta/lambda) = 3.10 cm at 30 keV
% absorption_length.value = 30.9;
% absorption_length.unit = 'cm';
%
% scan_time.value = 20 ;
% scan_time.unit = 's' ;
%
% object_length.value = 1.0 ;
% object_length.unit = 'mm';
%
% object_density.value  = 0.997 ;
% object_density.unit  = 'g / ccm';
% object_density.unit_equivalent = '1000 kg / m^3' ;
%
% container_length.value = 12 ; % mm
% container_length.unit = 'mm';
%
% %% P05 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% %% absorber
%
% % CVD 300 micron
% % Cu 50 micron
% % transmission at 30000 eV: 0.59
%
% % glassy / vitreous carbon / Glaskohlenstoff
%
% absorber(1).type = 'GC';
% % density
% absorber(1).density.value = [1.4 1.5];
% absorber(1).density.unit = 'g / ccm';
% % thickness
% absorber(1).thickness.value = 0.004;
% absorber(1).thickness.unit = 'm';
%
% % energy
% p05.energy.value = [14400 3000ingrid.fuchs@hzg.de0 ];
% p05.energy.unit = 'eV';
%
% % max flux at sample position
% p05.flux.value = [7e13, 525542981736 / 3 ];
% p05.flux.unit = 'photons / s';
%
% % spot size at sample position
% p05.area.value = [6.000 * 3.0, 1.0 * 1.0];
% p05.area.unit = 'mm^2';
%
% % flux density
% p05.flux_density.value = p05.flux.value ./ p05.area.value;
% p05.flux_density.unit = 'photons / s / mm^2';
%
% %% Materials
%
% % sample environment
%
% transmission.PEEK = 0.85; %extrapolated guess from 30
%
% % average bone+PEEK
% transmission.data.bone_PEEK = 0.61;
%
% %% Bone
%
% % formulas: https://physics.nist.gov/PhysRefData/XrayMassCoef/chap2.html
% % constants: https://physics.nist.gov/PhysRefData/XrayMassCoef/tab2.html
%
% bone.thickness.value = 0.001;
% bone.thickness.unit = 'm';
% bone.density.value = 1.920E+00 / 1000 / 0.01^3 ;
% bone.density.unit = 'kg / m^3';
% bone.density.unit_equiv = 'g / cm^3 * 1000 * 0.01^3';
%
% % Energy / MeV
% energy.value = [1.00000E-02 1.50000E-02 2.00000E-02 3.00000E-02 4.00000E-02 5.00000E-02 6.00000E-02 8.00000E-02];
% energy.unit = 'MeV';
%
% % X-ray mass attenuation coefficient mu/rho / (cm^2/g)
% mass_att_coeff.value = [ 2.851E+01 9.032E+00 4.001E+00 1.331E+00 6.655E-01 4.242E-01 3.148E-01 2.229E-01];
% mass_att_coeff.unit = 'cm^2 / g';
% mass_att_coeff.def = 'mu / rho';
%
% % transmission I/I_0 = exp( - mu * t )
% trans = exp( - mass_att_coeff.value * bone.density.value * bone.thickness.value * 100 );
% mu_over_rho_35000eV = mean( mass_att_coeff.value(4:5) );
% transmission_bone_35000eV = exp( - mu_over_rho_35000eV * bone.density.value * bone.thickness.value*100 );
% bone.volume.value = (bone.thickness.value)^3;
% bone.volume.unit = 'm^3';
% bone.mass.value = bone.density.value * bone.volume.value;
% bone.mass.unit = 'kg';
%
% absorption = 1 - transmission_bone_35000eV;
%
% % exposure
% exposure_time.value = [1200 * 0.15, 3000 * 0.1];
% exposure_time.unit =  's';
%
% energy.value = 34000;
% energy.unit = 'eV';
%
% %% Dose
% flux_density_30000eV = 1;
% dose = transmission.PEEK * absorption * flux_density_30000eV * (bone.thickness.value*1000)^2 * exposure_time.value * energy * eV_to_J / bone.mass.value;
%
% fprintf( '\ndose : %g\b kGy', dose / 1000 )
%
% % Attenuation
% out.containerTransmission = exp( - (container_length.value -object_length.value)/2 / absorption_length.value) ;
% out.objectAbsorption = 1-exp( - object_length.value/ absorption_length.value) ;
% % Number of photons
% out.absorbedPhotons = out.containerTransmission * out.objectAbsorption * fluxDensity_ph_per_s_mm2 * scan_time.value  * (object_length.value)^2;
% % Cubic mass of object
% out.cubeMass_g = object_length.value^3 * objectDensity_g_per_ccm / 1000 ;
% % Dose
% out.dose_Gy = out.absorbedPhotons * energy.value*10^3 * eV_to_J ...
%     / ( out.cubeMass_g / 1000 );

