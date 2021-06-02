% One gray is the absorption of one joule of energy, in the form of
% ionizing radiation, per kilogram of matter.
%
% 1 Gy = 1 J / kg = m^2 / s^2

%ca
%clear all
sample_type = 'real';
sample_type = 'cylinder';
num_proj_phan = 1200;
%scintillator_factor = cdwo300.absorption.value;

eV_to_J = 1 / 6.24e18;

bin = 4;
energy = 35000:5000:45000;
fprintf( ' energy: min : %g keV, max : %g keV, steps : %u', min( energy )/1000, max( energy )/1000, numel( energy ) )
voxel_size.value = bin * 2 * 3.115e-6;
voxel_size.unit = 'm';
fprintf( '\n voxel_size : %g micron', voxel_size.value * 1000^2)
phantom_diameter = 7e-3;
exp_time_per_image.value = 34e-3;
exp_time_per_image.unit = 's';

font_size = 18;
line_width = 6;
fig_path = '/asap3/petra3/gpfs/p05/2019/data/11006704/processed/dose/';
CheckAndMakePath( fig_path )

%% Scintillator
cdwo100 = cadmium_tungstate(energy,100e-6);
cdwo300 = cadmium_tungstate(energy,300e-6);
Y = [cdwo100.absorption.value; cdwo300.absorption.value];
if exist( 'h1' , 'var' ) && isvalid( h1 )
    figure(h1)
else
    h1 = figure( 'Name', 'Absorption Scintillator Cadmium Tungstate' );
end
p = plot( energy / 1000, Y );
legend( {'100 micron', '300 micron'}, 'FontSize',font_size )
title( 'scintillator efficiency' )
xlabel( 'energy / keV', 'FontSize',font_size )
ylabel( 'absorption', 'FontSize',font_size )
set( p ,'LineWidth', line_width )
ax = gca;
ax.FontSize = font_size; 
axis tight
saveas( h1, sprintf( '%s%s.png', fig_path, regexprep( h1.Name, '\ |:', '_') ) );

%% PEEK cylinder
peek05 = PEEK( energy, 5e-3 );
peek10 = PEEK( energy, 10e-3 );
%Y = [peek05.transmission.value; peek10.transmission.value];
% if exist( 'hpeek' , 'var' ) && isvalid( hpeek )
%     figure(hpeek)
% else
%     hpeek = figure( 'Name', 'Transmission PEEK scintillator' );
% end
% plot( energy / 1000, Y)
% legend( {'5 mm', '10 mm'} )
% xlabel( 'energy / keV' )

%% Bone sample
switch sample_type
    case 'real'
        scan_path = '/asap3/petra3/gpfs/p05/2017/data/11003440/processed/syn32_99R_Mg10Gd_4w/segmentation/2017_11003440_syn32_99R_Mg10Gd_4w_labels';
        if ~exist( 'vol', 'var' )
            fprintf( '\n' )
            vol = read_images_to_stack( scan_path );
            domain( vol(:), 1, 'original volume  ' )
            vol = single( vol );
            vol = vol( 1:bin:end, 1:bin:end, 1:bin:end );
%             vol = Binning( vol, bin ) / bin^3;
%             domain( vol(:), 1, 'single conversion' )
%             vol = permute( vol, [1 2 3] );
        end
        fprintf( '\nvolume shape : %u %u %u', size( vol ) )
        vol_bone = double( vol == 205 );
        
        %sc =  (vol == 105 | vol == 5);
    case 'cylinder'
        vol = zeros( [300 300], 'double' );
        d12 = phantom_diameter / voxel_size.value  / 2;
        x = (1:size(vol,1) ) - round( size(vol,1) / 2 );
        y = (1:size(vol,2) ) - round( size(vol,2) / 2 );
        [xx,yy] = meshgrid( x, y );
        m = sqrt( xx.^2 + yy.^2 );
        vol( m < d12 ) = 1;
        vol_bone  = repmat( vol, [1 1 250] );
        clear vol
end
bc = bone_cortical(energy,2e-3);
figure( 'Name', 'Sample' )
subplot( 1, 3, 1)
im = squeeze( vol_bone(:,:, round( size(vol_bone,3)/2)));
imwrite( im, sprintf( '%ssample_%s_noFig.png', fig_path, sample_type ) );
imsc( im ) 
axis equal tight
subplot( 1, 3, 2)
imsc( squeeze( vol_bone(:,round( size(vol_bone,2)/2),:)))
axis equal tight
subplot( 1, 3, 3)
imsc( squeeze( vol_bone(round( size(vol_bone,1)/2),:,:)))
axis equal tight

f = figure( 'Name', 'Sample hor cut' );
imsc( im )
title( ['phantom: ' sample_type], 'FontSize', font_size + 1  );
axis equal tight
set(gca,'xtick',[])
set(gca,'xticklabel',[])
set(gca,'ytick',[])
set(gca,'yticklabel',[])
saveas( f, sprintf( '%ssample_%s.png', fig_path, sample_type ) );

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
peek_cyl_rad_inner.value = 45e-3 / 2;
peek_cyl_rad_inner.unit = 'm';
peek_cyl_rad_outer.value = 55e-3 / 2;
peek_cyl_rad_outer.unit = 'm';
f = @(x) sqrt(peek_cyl_rad_outer.value^2 - x.^2) - sqrt(peek_cyl_rad_inner.value^2 - x.^2);
xx = size( vol_bone, 1 );
x = voxel_size.value * ( ( 1:xx ) - round( xx / 2 ) );
projected_thickness_peek.value = f( x );
projected_thickness_peek.unit = 'm';
fprintf( '\nPEEK cylinder: ' )
fprintf( '\n outer diameter : %g mm', peek_cyl_rad_outer.value * 1000 )
fprintf( '\n inner diameter : %g %s', peek_cyl_rad_inner.value * 1000)
fprintf( '\n projected thicknes: [min, max] = [%.2g %.2g] mm', min( projected_thickness_peek.value ) * 1000, max( projected_thickness_peek.value ) * 1000 )
% if exist( 'hpeek_thick' , 'var' ) && isvalid( hpeek_thick )
%     figure(hpeek_thick)
% else
%     hpeek_thick = figure( 'Name', 'PEEK cylinder thickness' );
% end
% plot( x * 1000, projected_thickness_peek.value * 1000)
% xlabel( 'radius / mm' )
% ylabel( 'thickness / mm' )

%% Transmission water
water_cyl_rad_inner.value = 10e-3 / 2;
water_cyl_rad_inner.unit = 'm';
water_cyl_rad_outer.value = 12e-3 / 2;
water_cyl_rad_outer.unit = 'm';
f = @(x) sqrt(water_cyl_rad_outer.value^2 - x.^2) - sqrt(water_cyl_rad_inner.value^2 - x.^2);
xx = size( vol_bone, 1 );
x = voxel_size.value * ( ( 1:xx ) - round( xx / 2 ) );
projected_thickness_water.value = f( x );
projected_thickness_water.unit = 'm';
fprintf( '\nWATER cylinder: ' )
fprintf( '\n outer diameter : %g mm', water_cyl_rad_outer.value * 1000 )
fprintf( '\n inner diameter : %g %s', water_cyl_rad_inner.value * 1000)
fprintf( '\n projected thicknes: [min, max] = [%.2g %.2g] mm', min( projected_thickness_water.value ) * 1000, max( projected_thickness_water.value ) * 1000 )

%% Dose %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
angles = (0:num_proj_phan-1)/num_proj_phan * pi;
projected_thickness_bone = permute( astra_make_sino_3D( vol_bone, angles ), [ 3 1 2] ) * voxel_size.value;
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
fd.value =  [5.4 5.4] * 1e10 / (0.0065 * 0.0048);
fd.unit = 'photons / s / m^2';
fd_energy.value = [30 45] * 1e3;
fd_energy.unit = 'eV';
flux_density.value = interp1( fd_energy.value, fd.value, energy, 'linear', 'extrap' );
flux_density.unit = 'photons / s / m^2';

area = size( projected_thickness_bone, 1) * size( projected_thickness_bone, 2) * voxel_size.value^2;
fprintf( '\ndetector area : %g mm x %g mm', size( projected_thickness_bone, 1) * voxel_size.value * 1000, size( projected_thickness_bone, 2) * voxel_size.value * 1000 )
%flux_density.value * area;
flux.value = 1 * 1e11 * ones( [1  numel( energy ) ] ) ;

figure( 'Name', sprintf( 'Flux' ) );
plot( energy / 1000, flux.value )
xlabel( 'energy / keV' )
ylabel( 'photons / s' )
drawnow

%flux.value = 7.2e10 * ones( numel( energy ) );
flux.unit = 'photons / s';

bone_density = bc.density.value;
mac_bone = bc.mass_att_coeff.value;

peek_density = peek05.density.value;
mac_peek = peek05.mass_att_coeff.value;

water_density = 1.0 * 1000; % water density in kg / m^3
[energy_wl, mac_wl] = read_nist_txt( 'water_liquid' );
mac_water = interp1( energy_wl, mac_wl, energy );

%% FIGURE skin entrance dose
entrance_doserate_water = eV_to_J * flux_density.value .* mac_water .* energy;
%entrance_doserate_water = eV_to_J * 2.1e11*1e6 .* mac_water .* energy;
entrance_dose_water = entrance_doserate_water * total_scan_time;
name = sprintf( 'skin entrance dose water');
f = figure( 'Name', name );
p = plot( energy / 1000, entrance_dose_water );
xlabel( 'energy / keV', 'FontSize', font_size )
ylabel( 'dose / Gy', 'FontSize', font_size )
set( p ,'LineWidth', line_width )
title( name , 'FontSize', font_size );
ax = gca;
ax.FontSize = font_size; 
axis tight

%% FIGURE skin entrance bone
entrance_doserate_bone = eV_to_J * flux_density.value .* mac_bone.* energy;
entrance_dose_bone = entrance_doserate_bone* total_scan_time;
name = sprintf( 'skin entrance dose bone');
f = figure( 'Name', name );
p = plot( energy / 1000, entrance_dose_bone );
xlabel( 'energy / keV', 'FontSize', font_size )
ylabel( 'dose / Gy', 'FontSize', font_size )
set( p ,'LineWidth', line_width )
title( name , 'FontSize', font_size );
ax = gca;
ax.FontSize = font_size; 
axis tight

%% Loop
fprintf( '\n Start loop over energy and projections' )
exp_time = exp_time_per_image.value;
fprintf( '\n energy step (%u): ', numel( energy ) )
for kk = 1:numel( energy )
    fprintf( ' %u', kk)
    mac_bone_kk = mac_bone(kk);
    mac_peek_kk = mac_peek(kk);
    mac_water_kk = mac_water(kk);
    energy_kk = energy(kk);
    t = exp( - peek_density * projected_thickness_peek.value * mac_peek_kk );
    %t = t .* exp( - water_density * projected_thickness_water.value * mac_water_kk );
    t = repmat( t, [300 1] );
    % Flux per pixel
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

% Correction factor of two
% flux measured with fast shuter is smaller than without, however
% measurement with flux only accounts for opening flux reduction but not
% for closing ( from 7.1 without to 5.4 with opening flank fast shutter, so
% approximately 3.7 with opening and closing flank) AND flux was measured
% with 0.1 s instead of 0.05 because of  technical problems (from 5.4 at
% 0.1 ms to 3.7 at 0.05 ms)
dose.value = absorbed_energy_per_scan / bone_mass / 2;
dose.unit = 'Gy';

%% FIGURE dose
f = figure( 'Name', sprintf( 'Dose %s: Cortical bone', sample_type) );
p = plot( energy / 1000, dose.value / 1000 );
xlabel( 'energy / keV', 'FontSize', font_size )
ylabel( 'dose / kGy', 'FontSize', font_size )
set( p ,'LineWidth', line_width )
title( sprintf('dose simulation: %s', sample_type ) , 'FontSize', font_size );
ax = gca;
ax.FontSize = font_size; 
axis tight
saveas( f, sprintf( '%s%s.png', fig_path, regexprep( f.Name, '\ |:', '_') ) );
switch sample_type
    case 'real'
        dose_real = dose;
    case 'cylinder'
        dose_phan = dose;
end

%% FIGURE efficiency scaled dose
f = figure( 'Name', sprintf( 'scintillator efficiency scaled dose %s: cortical bone', sample_type) );
x = energy / 1000;
[~,ind] = min( abs(energy - 45000) );
%fsc = @(x) x - min(x(:))  + 1;
scintillator_factor100 = ( cdwo100.transmission.value);
scintillator_factor300 = ( cdwo300.transmission.value);
% Normalize with scaling 1 at 45 keV where flux was measured
fsc = @(x) x - scintillator_factor300(ind)  + 1;
Y = dose.value / 1000 .* [fsc(scintillator_factor100); fsc(scintillator_factor300)];
p = plot( x, Y );
xlabel( 'energy / keV', 'FontSize', font_size )
ylabel( 'effective dose / kGy','FontSize', font_size )
legend( {'100 micron', '300 micron'}, 'FontSize', font_size, 'Location','north' )
set( p ,'LineWidth', line_width )
title( sprintf('effective dose simulation: %s', sample_type ) , 'FontSize', font_size );
ax = gca;
ax.FontSize = font_size; 
axis tight
saveas( f, sprintf( '%s%s.png', fig_path, regexprep( f.Name, '\ |:', '_') ) );

%% FIGURE dose real plus phan
if exist( 'dose_real', 'var' ) && exist( 'dose_phan', 'var')
    f = figure( 'Name', sprintf( 'Dose: Cortical bone phantom vs real'));
    x = energy / 1000;
    Y = [dose_phan.value / 1000; dose_real.value / 1000];
    p = plot( x, Y );
    xlabel( 'energy / keV', 'FontSize', font_size )
    ylabel( 'dose / kGy', 'FontSize', font_size )
    set( p ,'LineWidth', line_width )
    title( sprintf('dose simulation' ) , 'FontSize', font_size );
    legend( {'cylinder', 'real'}, 'FontSize', font_size )
    ax = gca;
    ax.FontSize = font_size;
    axis tight
    saveas( f, sprintf( '%s%s.png', fig_path, regexprep( f.Name, '\ |:', '_') ) );
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf( '\n' )
