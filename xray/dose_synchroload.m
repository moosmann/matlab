% One gray is the absorption of one joule of energy, in the form of
% ionizing radiation, per kilogram of matter.
%
% 1 Gy = 1 J / kg = m^2 / s^2

%sample_type = 'real';
sample_type = 'phantom';
%scintillator_factor = cdwo300.absorption.value;

eV_to_J = 1 / 6.24e18;

bin = 4;
energy = 30000:2000:50000;
fprintf( ' energy: min : %g keV, max : %g keV, steps : %u', min( energy )/1000, max( energy )/1000, numel( energy ) )
voxel_size.value = bin * 2 * 3.115e-6;
voxel_size.unit = 'm';
fprintf( '\n voxel_size : %g micron', voxel_size.value * 1000^2)

exp_time_per_image.value = 300e-3;
exp_time_per_image.unit = 's';

%% Scintillator

cdwo100 = cadmium_tungstate(energy,100e-6);
cdwo300 = cadmium_tungstate(energy,300e-6);
scintillator_factor = cdwo100.transmission.value;
Y = [cdwo100.absorption.value; cdwo300.absorption.value];
if exist( 'h1' , 'var' ) && isvalid( h1 )
    figure(h1)
else
    h1 = figure( 'Name', 'Absorption Scintillator' );
end
plot( energy / 1000, Y)
legend( {'100 micron', '300 micron'} )
xlabel( 'energy / keV' )

%% PEEK cylinder
peek05 = PEEK( energy, 5e-3 );
peek10 = PEEK( energy, 10e-3 );
Y = [peek05.transmission.value; peek10.transmission.value];
if exist( 'hpeek' , 'var' ) && isvalid( hpeek )
    figure(hpeek)
else
    hpeek = figure( 'Name', 'Transmission PEEK scintillator' );
end
plot( energy / 1000, Y)
legend( {'5 mm', '10 mm'} )
xlabel( 'energy / keV' )

%% Bone sample
switch sample_type
    case {1, 'real'}
        scan_path = '/asap3/petra3/gpfs/p05/2017/data/11003440/processed/syn32_99R_Mg10Gd_4w/segmentation/2017_11003440_syn32_99R_Mg10Gd_4w_labels';
        if ~exist( 'vol', 'var' )
            fprintf( '\n' )
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
    case {2, 'phantom'}
        vol = zeros( [300 300], 'single' );
        d12 = 5e-3 / voxel_size.value  / 2;
        x = (1:size(vol,1) ) - round( size(vol,1) / 2 );
        y = (1:size(vol,2) ) - round( size(vol,2) / 2 );
        [xx,yy] = meshgrid( x, y );
        m = sqrt( xx.^2 + yy.^2 );
        vol( m < d12 ) = 1;
        vol_bone  = repmat( vol, [1 1 250] );
        clear vol
end
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

%% Transmission PEEK cylinder
peek_cyl_diam_inner.value = 45e-3 / 2;
peek_cyl_diam_inner.unit = 'm';
peek_cyl_diam_outer.value = 50e-3 / 2;
peek_cyl_diam_outer.unit = 'm';
f = @(x) sqrt(peek_cyl_diam_outer.value^2 - x.^2) - sqrt(peek_cyl_diam_inner.value^2 - x.^2);
xx = size( vol_bone, 1 );
x = voxel_size.value * ( ( 1:xx ) - round( xx / 2 ) );
projected_thickness_peek.value = f( x );
projected_thickness_peek.unit = 'm';
fprintf( '\nPEEK cylinder: ' )
fprintf( '\n outer diameter : %g mm', peek_cyl_diam_outer.value * 1000 )
fprintf( '\n inner diameter : %g %s', peek_cyl_diam_inner.value * 1000)
fprintf( '\n projected thicknes: [min, max] = [%.2g %.2g] mm', min( projected_thickness_peek.value ) * 1000, max( projected_thickness_peek.value ) * 1000 )
if exist( 'hpeek_thick' , 'var' ) && isvalid( hpeek_thick )
    figure(hpeek_thick)
else
    hpeek_thick = figure( 'Name', 'PEEK cylinder thickness' );
end
plot( x * 1000, projected_thickness_peek.value * 1000)
xlabel( 'radius / mm' )
ylabel( 'thickness / mm' )

%% Dose
% Projected bone thickness
projected_thickness_bone = permute( astra_make_sino_3D( vol_bone ), [ 3 1 2] ) * voxel_size.value;
fprintf( '\n projections shape : %u %u %u', size( projected_thickness_bone ) )
num_proj = size( projected_thickness_bone, 3 );
fprintf( '\n bone projected thickness: max = %g mm, avg = %g mm', max(projected_thickness_bone(:))*1000, mean( projected_thickness_bone(:)) *1000)

% Scan time
total_scan_time = num_proj * exp_time_per_image.value;
fprintf( '\n total scan time : %g s', total_scan_time )

% Preallocation
absorbed_energy_per_image = zeros( [1, num_proj] );
absorbed_energy_per_scan = zeros( [1, numel( energy )] );

% Flux values from previous experiment
fd.value = [2.5 3.7 5.2] *1e17;
fd.unit = 'photons / s / m^2';
fd_energy.value = [30 40 45] * 1e3;
fd_energy.unit = 'eV';
flux_density.value = interp1( fd_energy.value, fd.value, energy, 'linear', 'extrap' );
flux_density.unit = 'photons / s / m^2';

area = size( projected_thickness_bone, 1) * size( projected_thickness_bone, 2) * voxel_size.value^2;
fprintf( '\ndetector area : %g mm x %g mm', size( projected_thickness_bone, 1) * voxel_size.value * 1000, size( projected_thickness_bone, 2) * voxel_size.value * 1000 )
flux.value = flux_density.value * area;

figure( 'Name', sprintf( 'Flux' ) );
plot( energy / 1000, flux.value )
xlabel( 'energy / keV' )
ylabel( 'photons / s' )
drawnow

%flux.value = 7.2e10 * ones( numel( energy ) );
flux.unit = 'photons / s';

fprintf( '\n Start loop over energy and projections' )
bone_density = bc.density.value;
mac_bone = bc.mass_att_coeff.value;
peek_density = peek05.density.value;
mac_peek = peek05.mass_att_coeff.value;
exp_time = exp_time_per_image.value;
fprintf( '\n energy step (%u): ', numel( energy ) )
for kk = 1:numel( energy )
    fprintf( ' %u', kk)
    mac_bone_kk = mac_bone(kk);
    mac_peer_kk = mac_peek(kk);
    energy_kk = energy(kk);
    t = exp( - peek_density * projected_thickness_peek.value * mac_peer_kk );
    t = repmat( t, [300 1] );
    f_kk = flux.value(kk) / size( projected_thickness_bone, 1 ) / size( projected_thickness_bone, 2 );
    for ll = 1:num_proj
        % absorption image
        a = 1 - exp( - bone_density * projected_thickness_bone(:,:,ll) * mac_bone_kk );
        % absorbed energy of full image
        absorbed_energy_per_image(ll) = sum( sum( a * t * f_kk * exp_time .* energy_kk * eV_to_J ) );
    end
    % absorbed energy of full tomo scan summing over all images
    absorbed_energy_per_scan(kk) = sum( absorbed_energy_per_image );
end
fprintf( '\n Loop finished' )

dose.value = absorbed_energy_per_scan / bone_mass;
dose.unit = 'Gy';

figure( 'Name', sprintf( 'Dose %s: Cortical bone', sample_type) );
plot( energy / 1000, dose.value / 1000 )
xlabel( 'energy / keV' )
ylabel( 'dose / kGy' )

%% efficiency scaled dose
figure( 'Name', sprintf( 'Dose*scintillator efficiency %s: Cortical bone', sample_type) );
x = energy / 1000;

scintillator_factor = scintillator_factor - min( scintillator_factor(:) ) + 1;
y = dose.value / 1000 .* scintillator_factor;
plot( x, y )
xlabel( 'energy / keV' )
ylabel( 'dose / kGy' )


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

