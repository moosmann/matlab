

clear all
tic;
% DLS data
%p1 = '/asap3/petra3/gpfs/external/2019/data/50000258/processed/resampled/113769';
p1 = '/asap3/petra3/gpfs/external/2019/data/50000258/processed/resampled/113734';
%p2 = '/asap3/petra3/gpfs/external/2019/data/50000258/processed/resampled/113734';
e1 = 28e3; %eV, pink beam diamond;

% PETRA III data
mpath = '/asap3/petra3/gpfs/p05/2016/data/11001978/processed/segmentation/matlab/scans_resampled.mat';
load( mpath );
%p1 = '/asap3/petra3/gpfs/p05/2016/data/11001978/processed/segmentation/scans/cpd_Mg10Gd_E33998eV_cam3056x3056_pix2420nm77R_2018_11004263_syn019_77R_Mg10Gd_8w';
%p1 = '/asap3/petra3/gpfs/p05/2016/data/11001978/processed/segmentation/scans/cpd_Mg5Gd_E34014eV_cam3056x1621_pix2427nm77L_2017_11003950_syn22_77L_Mg5Gd_8w';
%p1 = '/asap3/petra3/gpfs/p05/2016/data/11001978/processed/segmentation/scans/cpd_Mg5Gd_E34014eV_cam3056x1621_pix2427nm80L_2017_11003950_syn22_80L_Mg5Gd_8w';
p2 = '/asap3/petra3/gpfs/p05/2016/data/11001978/processed/segmentation/scans/cpd_Mg5Gd_E34014eV_cam3056x1621_pix2427nm88R_2017_11003950_syn22_88R_Mg5Gd_4w';
p2 = '/asap3/petra3/gpfs/p05/2016/data/11001978/processed/segmentation/scans/cpd_Mg10Gd_E38374eV_cam3056x3056_pix2403nm_82R_2017_11003440_syn21_82R_Mg10Gd_8w';

[~, name] = fileparts( p2 );
m = contains( {scans_resampled.name}, name); 
scan = scans_resampled(m);
if isempty( scan.scan_structs )
    e2 = scan.scans.energy;
    mat2 = scan.scans.material;
else
    e2 = scan.scan_structs(1).energy;
    mat2 = scan.scan_structs(1).material;
end


% Energies
erange = [e1 e2];
fprintf( '\n energies / keV: ' )
fprintf( ' %.1f', erange / 1000 )

% Bone
[bc, mac_bc(1)] = bone_cortical( e1 );
[~, mac_bc(2)] = bone_cortical( e2 );
rho_bc = bc.density.value;
ac_bc = rho_bc * mac_bc;
fprintf( '\nCortical bone:')
fprintf( '\n density / (%s) : %g', bc.density.unit, rho_bc );
fprintf('\n mass attenuation coefficient / (%s) :',bc.mass_att_coeff.unit )
fprintf( ' %g', mac_bc )
fprintf('\n attenuation coefficient / m :' )
fprintf( ' %g', ac_bc * 5e-6)

% Mg-Gd alloy
% Mg5Gd
[mg5gd, mac_mg5gd(1)] = magnesium_gadolinium( e1, 5 );
[~, mac_mg5gd(2)] = magnesium_gadolinium( e2, 5 );
rho_mg5gd = mg5gd.density.value;
ac_mg5gd = rho_mg5gd * mac_mg5gd;
fprintf( '\nMg5Gd:' )
fprintf( '\n density / (%s) : %g', mg5gd.density.unit, rho_mg5gd );
fprintf('\n mass attenuation coefficient / (%s) :',mg5gd.mass_att_coeff.unit )
fprintf( ' %g', mac_mg5gd )
fprintf('\n attenuation coefficient / m :' )
fprintf( ' %g', ac_mg5gd * 5e-6)
% Mg10Gd
[mg10gd, mac_mg10gd(1)] = magnesium_gadolinium( e1, 10 );
[~, mac_mg10gd(2)] = magnesium_gadolinium( e2, 10 );
rho_mg10gd = mg10gd.density.value;
ac_mg10gd = rho_mg10gd * mac_mg10gd;
fprintf( '\nMg10Gd:' )
fprintf( '\n density / (%s) : %g', mg10gd.density.unit, rho_mg10gd );
fprintf('\n mass attenuation coefficient / (%s) :',mg10gd.mass_att_coeff.unit )
fprintf( ' %g', mac_mg10gd )
fprintf('\n attenuation coefficient / m :' )
fprintf( ' %g', ac_mg10gd * 5e-6)

fprintf( '\n' )
%

d1 = dir( [p1 '/*tif' ] );
d2 = dir( [p2 '/*tif' ] );

ni1 = numel( d1 );
ni2 = numel( d2 );

n1 = round( 0.5 * ni1  );
n2 = round( 0.5 * ni2 );

filename = [d1(n1).folder filesep d1(n1).name];
im1 = imread( filename );

filename = [d2(n2).folder filesep d2(n2).name];
im2 = imread( filename );

% Compression
parameter = [0.01 0.01];
verbose = 0;
[tlow1, thigh1] = compression( im1, 'histo', parameter, verbose );
%im1 = (im1 - tlow1) ./ (thigh1 - tlow1);
domain( im1 )
im1( im1 < tlow1 ) = tlow1;
im1( im1 > thigh1 ) = thigh1;
domain( im1 );
[tlow2, thigh2] = compression( im2, 'histo', parameter, verbose );
%im2 = (im2 - tlow2) ./ (thigh2 - tlow2);
im2( im2 < tlow2 ) = tlow2;
im2( im2 > thigh2 ) = thigh2;

% Plot images
imr = @(im) sprintf( 'min/max: %f/%f\nmean/std:%f/%f', min(im(:)), max(im(:)), mean(im(:)), std(im(:)) );

subplot( 2, 2, 1 )
imsc( im1 )
title( sprintf( 'DLS\n%s', imr(im1) ) )
axis equal tight

subplot( 2, 2, 2 )
imsc( im2 )
title( sprintf( 'PETRA III\n%s', imr(im2) ) )
axis equal tight

% Plot histograms
subplot( 2, 2, 3 )
num_bins = 1024;
bl1 = min( [ min(im1(:)), min(im2(:)) ] );
bl2 = max( [ max(im1(:)), max(im2(:)) ] );
bin_limits = [bl1 bl2];
[N1, edges1] = histcounts( im1, num_bins, 'BinLimits', bin_limits );
[N2, edges2] = histcounts( im2, num_bins, 'BinLimits', bin_limits );
X1 = ( edges1(2:end) + edges1(1:end-1) ) / 2;
Y1 = [N1; N2];
plot( X1, Y1 )
axis tight
lab_pet = sprintf( 'PETRA III %s', mat2 );
legend( { 'DLS', lab_pet }, 'Location','NorthWest' )
title( 'Same limits' )

%%
subplot( 2, 2, 4 )
[N1, edges1] = histcounts( im1, num_bins );
[N2, edges2] = histcounts( im2, num_bins );
X2 = 1:num_bins;
Y2 = [N1; N2]';
plot( X2, Y2 )
axis tight
lab_pet = sprintf( 'PETRA III %s', mat2 );
legend( { 'DLS', lab_pet }, 'Location','NorthWest' )
title( 'Individual limits' )


figure( 'Name', 'Histogram: data + tabulated values' )
%% 
num_bins = 1024;
bl1 = min( [ min(im1(:)), min(im2(:)) ] );
bl2 = max( [ max(im1(:)), max(im2(:)) ] );
bin_limits = [bl1 bl2];
[N1, edges1] = histcounts( im1, num_bins, 'BinLimits', bin_limits );
[N2, edges2] = histcounts( im2, num_bins, 'BinLimits', bin_limits );
X1 = ( edges1(2:end) + edges1(1:end-1) ) / 2;
Y1 = [N1; N2];
plot( X1, Y1 )
title( 'Same limits' )
e1l = sprintf( ' %.0fkeV', e1/1000);
e2l = sprintf( ' %.0fkeV', e2/1000);
scf = 3.5e-6;
lha = 'left';
xline( scf * ac_bc(1),'-', ['bone' e1l], 'LabelHorizontalAlignment', lha, 'LabelVerticalAlignment', 'top');
xline( scf * ac_bc(2),'-', ['bone' e2l], 'LabelHorizontalAlignment', lha, 'LabelVerticalAlignment', 'top');
xline( scf * ac_mg5gd(1),'-', ['Mg5Gd' e1l], 'LabelHorizontalAlignment', lha, 'LabelVerticalAlignment', 'middle' );
xline( scf * ac_mg5gd(2),'-', ['Mg5Gd' e2l], 'LabelHorizontalAlignment', lha, 'LabelVerticalAlignment', 'middle' );
xline( scf * ac_mg10gd(1),'-', ['Mg10Gd' e1l], 'LabelHorizontalAlignment', lha, 'LabelVerticalAlignment', 'top' );
xline( scf * ac_mg10gd(2),'-', ['Mg10Gd' e2l], 'LabelHorizontalAlignment', lha, 'LabelVerticalAlignment', 'top' );
lab_pet = sprintf( 'PETRA III %s', mat2 );
legend( { 'DLS', lab_pet }, 'Location','NorthWest' )
axis tight

if 0
    % Volumes
    v1 = read_images_to_stack( p1, 1, '*tif', [], 1, 1 );
    v2 = read_images_to_stack( p2, 1, '*tif', [], 1, 1 );
    
    %% Histo
    hc1 = histcounts( v1(:), num_bins );
    hc2 = histcounts( v2(:), num_bins );
    
    %% Plot histo
    figure('Name', 'Volume histogram' )
    X2 = 1:num_bins;
    Y2 = [normat( hc1 ); normat( hc2 )]';
    plot( X2, Y2 )
    axis tight
    legend( { 'DLS', 'PETRA III' } )    
    title( 'Individual limits' )
end