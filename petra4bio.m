%clear all

% emittance = eps_i = RMS-width of source * divergence eps_i = sigma_i *
% sigma_theta_i, RMS(f(x)) = sqrt( int_x1_x2 f(x)^2 dx )

% Brillicance = (spectral photon flux at 0.1% bandwidth at E) / (4 pi *
% eps_x * epx_y) = [ photons / s / rad^2 / m^2 / 0.1% bandwidth]

% Refractive index
delta_H20_30keV = 2.56114106E-07;
beta_H20_30keV = 1.06431745E-10;

%% Source size
% PETRA III
p3_sigma_h = 140e-6; % m
p3_sigma_v = 6e-6; % m
% PETRA IV
p4_sigma_h = 6e-6; % m
p4_sigma_v = p4_sigma_h; %m
% ESRF
esrf_sigma_h = 59e-6; % m. 51 FWHM: 90 micron
esrf_sigma_v = 8.3e-6; %m. 8.6 FWHM: 20 micron

%% emittance
p3_emitt_h = 1200e-12; % m * rad
p3_emitt_v = 10e-12; % m * rad
p4_emitt_h = 10e-12; % m * rad
p4_emitt_v = 10e-12; % m * rad

%% divergence
p3_div_h = p3_emitt_h / p3_sigma_h;
p3_div_v = p3_emitt_v / p3_sigma_v;
p4_div_h = p4_emitt_h / p4_sigma_h;
p4_div_v = p4_emitt_v / p4_sigma_v;

%% relative coherence
p3_coh_rel = 0.01; % per cent
p4_coh_rel = 0.25; % per cent

%% maximum brillicance
p3_brill_max = 1e20;
p4_brill_max = 1e22;

%% Energy
E = 30e3; % eV
dE = E * 1e-4;
lambda = E_to_lambda( E );
dlambda = lambda * dE / E;
w = E / reducedPlanckConstant;
dw = dE / reducedPlanckConstant;
f = E / PlanckConstant;
df = dE / PlanckConstant;
c = speedOfLight; %m/s

%% Distances
% distance source sample
p05_dist_source_dcm = 50.9; % m
p05_dist_source_sample = 82.7; % m
p05_dist_dcm_sample = p05_dist_source_sample - p05_dist_source_dcm; % m
id19_dist_source_sample = 145; % m

% longitudinal coherence lenght in vacuum v = c/n with n = 1
coh_l = speedOfLight * 2 * pi / dw; % m
coh_t_induced_by_coh_l = sqrt( 2 * p05_dist_source_sample * coh_l ); % m

coh.long__micron = coh_l * 1e6;
coh.transByLong__micron = coh_t_induced_by_coh_l * 1e6;

%% For printing
out.E__keV = E / 1000;
out.dE__keV = dE / 1000;
out.lambda__m = lambda * 1e10;
out.dlambda_m = dlambda * 1e10;
out.w__Hz = w;
out.dw__Hz = dw ;
out.f__Hz = f;
out.df__Hz = df;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N = 500;
num_angles = round( 1.6 * N );
sf = 1;
pixelsize = 1.4e-6;
dist_sample_detector = 0.5;
voxelsize = sf * pixelsize;
M = ceil( sqrt(2) * sf * N ) + 1;
gauss_noise_mean = 0.05;
delta_phase = 0.1;
coherence = 'p4';'esrf';'full';'p3';

phase_method = 'tie';
reg_par = 7;
bin_filt = 0.1;
cutoff_frequ = 2*pi;
phase_padding = 1;

vs = N / 2 * voxelsize;
min_x = - vs;
max_x = vs;
min_y = - vs;
max_y = vs;
min_z = - vs;
max_z = vs;
%num_angles = round( 2 * pi /2 * N);

angle_range = pi;
angles = angle_range * (0:num_angles - 1) / num_angles;
det_spacing_x = pixelsize;
det_spacing_y = pixelsize;
det_row_count = M;
det_col_count = M;

% Create phantom
phan = phantom3d( 'Modified Shepp-Logan' , N);
phan = FilterBlur( phan, 0, 1);
%phan = delta_H20_30keV * phan;
%phan = 0.25 * phan;
phan = delta_phase * phan / k0 / pixelsize / 0.25 / N;

% Forward projection
fprintf( '\nForward projection' );
k0 = 2 * pi / lambda;
vol_geom = astra_create_vol_geom( N, N, N, min_x, max_x, min_y, max_y, min_z, max_z );
proj_geom = astra_create_proj_geom('parallel3d', det_spacing_x, det_spacing_y, det_row_count, det_col_count, angles);
[sino_id, sino] = astra_create_sino3d_cuda( phan, proj_geom, vol_geom);
sino = k0 * sino;

% Forward propagation
fprintf( '\nForward propagation ' );
edp = [E, dist_sample_detector, pixelsize];
prop_padding = 1;
prop_method = 'symmetric';
sino_prop = zeros( size( sino ), 'like', sino );
% FWHM = 2 * sqrt( 2 * log(2) ) * std, choose var that FWHM is mean / 2
gauss_noise_var = (gauss_noise_mean / (4 * sqrt( 2 * log(2) ) ) )^2;
add_noise = @(im) im + imnoise( zeros(size( im )), 'gaussian', gauss_noise_mean, gauss_noise_var);
fshape = (1 + prop_padding) * [M M];
% Gaussian source convolution (blurring): VCZ in the far for an extended
% quasimonochromatic, incoherent source with Gaussian profile assuming weak
% phase variations and small absolute absorption
switch coherence
    case 'full'
        filt = 1;
    case 'p3'
        source_sigma_h_v = [p3_sigma_h, p3_sigma_v];
        dist_source_sample = p05_dist_source_sample;
    case 'p4'
        source_sigma_h_v = [p4_sigma_h, p4_sigma_v];
        dist_source_sample = p05_dist_source_sample;
    case 'esrf'
        source_sigma_h_v = [esrf_sigma_h, esrf_sigma_v];
        dist_source_sample = id19_dist_source_sample;
end
if ~strcmp( coherence, 'full' )
    filt = FilterGaussianSource(fshape, sigma_to_FWHM( source_sigma_h_v ), dist_source_sample, dist_sample_detector, pixelsize);
end
parfor nn = 1:size( sino, 2 )
    im = squeeze( sino(:,nn,:) );
    int = propagate( -im, 0*im, edp, prop_padding, prop_method, filt, 0);
    % Noise
    int = add_noise( int );
    sino_prop(:,nn,:) = int;
end

% Phase retrieval
sino_retr = zeros( size( sino_prop ), 'like', sino_prop );
im_shape = [ size( sino_prop, 1 ), size( sino_prop, 1 )];
[phase_filter, pha_appendix] = PhaseFilter( phase_method, (1 + phase_padding) * im_shape, edp, reg_par, bin_filt, cutoff_frequ, 'double');
parfor nn = 1:size( sino, 2 )
    im = padarray( squeeze( sino_prop(:,nn,:) ), phase_padding * im_shape, 'symmetric', 'post' );
    pha = -real( ifft2( phase_filter .* fft2( im ) ) );
    sino_retr(:,nn,:) = pha(1:im_shape(1), 1:im_shape(2));
end
%sino_retr = sino_retr + 0.5 / 10^-reg_par;

% Remove negative values
%sino_retr( sino_retr < 0 ) = 0;

% Ram-Lak Filter
fprintf( '\nFilter sino' );
fbp_filter_type = 'Ram-Lak';
fbp_filter_padding = 1;
fbp_filter_padding_method = 'symmetric';
fpb_filter_freq_cutoff = 1;
[N1, N2, N3] = size( sino_retr);
fbp_filt = iradonDesignFilter(fbp_filter_type, (1 + fbp_filter_padding) * N1, fpb_filter_freq_cutoff);
sino_pha_ret_fbpfilt = zeros( [N1, N2, N3], 'like', sino_retr );
parfor nn = 1:N2
    im = sino_retr(:,nn,:);
    im = padarray( im, fbp_filter_padding * [N1 0 0], fbp_filter_padding_method, 'post' );
    im = real( ifft( bsxfun(@times, fft( im, [], 1), fbp_filt), [], 1, 'symmetric') );
    sino_pha_ret_fbpfilt(:,nn,:) = im(1:N1,:,:);
end

% Backprojection
fprintf( '\nBack-projection' );
[reco_id, reco] = astra_create_backprojection3d_cuda( sino_pha_ret_fbpfilt, proj_geom, vol_geom);
% Correct inconsistent scaling in ASTRA
reco = 1 / k0 * (angles(end) - angles(1)) / 2 / num_angles / pixelsize / voxelsize^3 * reco;

% Error
diff_reco_phan = abs( reco - phan ) / (max(phan(:)) - min(phan(:)));

%coh.blur_sigma = blur_sigma;
%coh.blur_cutoff_frequ__1_m = blur_cutoff_frequ;
%coh.blur_cutoff_frequ_relativ = blur_cutoff_frequ * 2 * pixelsize;

sim.N = N;
sim.M = M;
sim.num_angles = num_angles;
sim.angle_first_last = [angles(1), angles(end)];

%% Print %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Parameters
disp( out )
disp( coh )
disp( sim )

% Ranges
domain( filt(:), 1, 'Fourier space blurring filter' );

domain( sino_prop(:), 1, 'intensities / propagated projections' )
domain( sino_pha_ret_fbpfilt(:), 1, 'Ram-Lak filtered retrieved phase maps' )

domain( sino(:), 1,      'projection: phase maps' )
domain( sino_retr(:), 1, 'retrieved phase maps' )

domain( diff_reco_phan(:), 1, 'diff: |reco - phan|' )
domain( phan(:), 1,           'Shepp-Logan phantom' )
domain( reco(:), 1,           'reconstruction     ' )

aclear;

%% Show %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nn = 1;

set( figure(nn),  'Name',  'Fourier space blurring filter' );
imsc( fftshift(filt));
axis tight equal; colorbar; nn = nn + 1;

set( figure(nn),  'Name',  'phantom: Shepp-Logan' );
sli = round( size( phan, 3) / 2 );
rphan = [ min(phan(:)) max(phan(:))];
imsc( phan(:,:,sli),  rphan);
axis tight equal; colorbar; nn = nn + 1;
title( sprintf( 'phan(:,:,%u)', sli ) )

set( figure(nn),  'Name',  'projection: phase map' );
imsc( squeeze(sino(:,1,:)), [0 max( sino(:))]);
axis tight equal; colorbar; nn = nn + 1;
title( sprintf( 'sino(:,1,:)' ) )

set( figure(nn),  'Name',  'projection: sino' );
sli = round( size( sino, 3) / 2 );
imsc( squeeze(sino(:,:,sli)), [0 max( sino(:))]);
axis tight equal; colorbar; nn = nn + 1;
title( sprintf( 'sino(:,:,%u)', sli ) )

set( figure(nn),  'Name',  'propagated projections: intensities' );
imsc( squeeze(sino_prop(:,1,:)) );
axis tight equal; colorbar; nn = nn + 1;
title( sprintf( '[sino prop](:,1,:)' ) )

set( figure(nn),  'Name',  'FT of propagated projections' );
im = squeeze(sino_prop(:,1,:));
im = log( 1 + abs( fftshift( fft2( SubtractMean( im ) ) ) ) );
imsc( im );
axis tight equal; colorbar; nn = nn + 1;
title( sprintf( 'FT[sino prop](:,1,:)' ) )


set( figure(nn),  'Name',  'retrieved phase maps' );
imsc( squeeze(sino_retr(:,1,:)) );
axis tight equal; colorbar; nn = nn + 1;
title( sprintf( '[sino retr](:,1,:)' ) )

set( figure(nn),  'Name',  'reconstruction' );
sli = round( size( reco, 3) / 2 );
imsc( squeeze(reco(:,:,sli)), rphan );
%imsc( squeeze(reco(:,:,sli)) );
axis tight equal; colorbar; nn = nn + 1;
title( sprintf( 'reco(:,:,%u)', sli ) )

set( figure(nn),  'Name',  'reconstruction: difference maps |reco - phan|' );
sli = round( size( reco, 3) / 2 );
imsc( squeeze(diff_reco_phan(:,:,sli)));
axis tight equal; colorbar; nn = nn + 1;
title( sprintf( '1 / max(phan(:)) * |reco - phan|(:,:,%u)', sli ) )
