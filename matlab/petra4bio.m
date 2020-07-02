%clear all
warning( 'off', 'MATLAB:imagesci:rtifc:missingPhotometricTag');

% emittance = eps_i = RMS-width of source * divergence eps_i = sigma_i *
% sigma_theta_i, RMS(f(x)) = sqrt( int_x1_x2 f(x)^2 dx )

% Brilliance = (spectral photon flux at 0.1% bandwidth at E) / (4 pi *
% eps_x * epx_y) = [ photons / s / rad^2 / m^2 / 0.1% bandwidth]

% Refractive index
delta_H20_30keV = 2.56114106E-07;
beta_H20_30keV = 1.06431745E-10;

%% Source size
% PETRA III
p3_sigma_h = 38e-6; % m
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
k0 = 2 * pi / lambda;

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

N = 1024; 
num_angles = round( 1.6 * N ); 
sf = 1;
pixelsize = 1e-6; 
dist_sample_detector = 0.5; 
voxelsize = sf * pixelsize;
M = ceil( sqrt(2) * sf * N ) + 1;
delta_phase = 0.5; 
gauss_noise_mean = 0.04; 
coherence = 'p4';'p3';'esrf';'full';

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
vol_geom = astra_create_vol_geom( N, N, N, min_x, max_x, min_y, max_y, min_z, max_z );
proj_geom = astra_create_proj_geom('parallel3d', det_spacing_x, det_spacing_y, det_row_count, det_col_count, angles);
[sino_id, sino] = astra_create_sino3d_cuda( phan, proj_geom, vol_geom);
sino = k0 * sino;

% Forward propagation
fprintf( '\nForward propagation ' );
edp = [E, dist_sample_detector, pixelsize];
prop_padding = 1;
prop_method = 'symmetric';
int = zeros( size( sino ), 'like', sino );
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
    im = propagate( -im, 0*im, edp, prop_padding, prop_method, filt, 0);    
    int(:,nn,:) = im;
end

% Noise
fprintf( '\nNoise' );

int_noise = zeros( size( int ), 'like', int );
parfor nn = 1:size( int, 2 )
    im = squeeze( int(:,nn,:) );    
    int_noise(:,nn,:) = add_noise( im );
end

% Phase retrieval
fprintf( '\nPhase retrieval' );
sino_retr = zeros( size( int_noise ), 'like', int_noise );
im_shape = [ size( int_noise, 1 ), size( int_noise, 1 )];
[phase_filter, pha_appendix] = PhaseFilter( phase_method, (1 + phase_padding) * im_shape, edp, reg_par, bin_filt, cutoff_frequ, 'double');
parfor nn = 1:size( sino, 2 )
    im = padarray( squeeze( int_noise(:,nn,:) ), phase_padding * im_shape, 'symmetric', 'post' );
    pha = -real( ifft2( phase_filter .* fft2( im ) ) );
    sino_retr(:,nn,:) = pha(1:im_shape(1), 1:im_shape(2));
end

% Ram-Lak Filter
fprintf( '\nFBP filter sino' );
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
diff_reco_phan = abs( reco - phan );
%coh.blur_sigma = blur_sigma;
%coh.blur_cutoff_frequ__1_m = blur_cutoff_frequ;
%coh.blur_cutoff_frequ_relativ = blur_cutoff_frequ * 2 * pixelsize;

sim.N = N;
sim.M = M;
sim.num_angles = num_angles;
sim.angle_first_last = [angles(1), angles(end)];

par.N = N;
par.num_angles = num_angles;
par.delta_phase = delta_phase;
par.gauss_noise_mean = 0.02;
par.dist_sample_detector = dist_sample_detector;
par.pixelsize = pixelsize;
par.coherence = coherence;
par.phase_method = phase_method;

%% In vivo data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dpath = '/asap3/petra3/gpfs/p05/2016/data/11001994/raw/szeb_74_00/';
d = imread([dpath 'dark_0000.tif']);
f = imread([dpath 'ref_0000.tif']);
p =imread([dpath 'proj_0000.tif']);
roi1 = round(size(d),1)/4:3*round(size(d),1)/4;
roi2 = round(size(d),2)/4:3*round(size(d),2)/4;
thresh = [0.05 0.02];
d = Binning( FilterPixel( d(roi1,roi2), thresh ) ) ;
p = Binning( FilterPixel( p(roi1,roi2), thresh ) ) ;
p = p -d;
f = Binning( FilterPixel( f(roi1,roi2), thresh ) ) ;
f = f-d;
int1 = p ./ f;
dm = mean2( d );
ds = std2( d );
pm = mean2( p );
ps = std2( p );
fm = mean2( f );
fs = std2( f );
intm = mean2( int1 );
%ishow(d),ishow(f),ishow(p),ishow(ii)

% Assume Poisson noise then mean(I) = var(I), but unknown conversion factor
% conversion factor given by mean(I) / var(I)
cnt_conv_fac = fm / fs^2;
cnt = fm * cnt_conv_fac;
% signal ( mean(I) = lambda ) to noise ( std(I)=sqrt(lambda) )
noise = sqrt( cnt );
noise_rel = noise / cnt;

%% Print %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Parameters
disp( out )
disp( coh )
disp( sim )

disp( par )

% Ranges
domain( filt(:), 0, 'Fourier space blurring filter' );

domain( int(:), 0, 'intensities without noise' )
domain( int_noise(:), 0, 'intensities with noise' )
domain( sino_pha_ret_fbpfilt(:), 0, 'Ram-Lak filtered retrieved phase maps' )

domain( sino(:), 0,      'projection: phase maps' )
domain( sino_retr(:), 0, 'retrieved phase maps' )

domain( diff_reco_phan(:)/max(phan(:)), 0, 'diff: |reco - phan|/max(phan)' )
domain( phan(:), 0,           'Shepp-Logan phantom' )
domain( reco(:), 0,           'reconstruction     ' )

aclear;

%% Show %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

parpath = sprintf('/asap3/petra3/gpfs/common/p05/jm/petra4bio/N%04u_sigHor%u/', N, p3_sigma_h*1e6 );
outpath = sprintf('%s%s/', parpath, coherence );
CheckAndMakePath( outpath );
save_png_par = @(im, name) imwrite( uint8(2^8*rescale(im)), [parpath name '.png'] );
save_png = @(im, name) imwrite( uint8(2^8*rescale(im)), [outpath name '.png'] );

nn = 1;

% Phantom
set( figure(nn),  'Name',  'phantom: Shepp-Logan' );
sli = round( size( phan, 3) / 2 );
im = phan(:,:,sli);
imsc( im );
axis tight equal; colorbar; nn = nn + 1;
title( sprintf( 'phan(:,:,%u)', sli ) )
name = sprintf( 'phantom_z%u', sli);
save_png( im, name)
% Tomo
set( figure(nn),  'Name',  'reconstruction' );
sli = round( size( reco, 3) / 2 );
im = squeeze( reco(:,:,sli) );
imsc( im );
axis tight equal; colorbar; nn = nn + 1;
title( sprintf( 'reco(:,:,%u)', sli ) )
name = sprintf( 'reco_tomo_z%u', 1);
save_png( im, name)

% simulated phase map
set( figure(nn),  'Name',  'projection: phase map' );
im = squeeze(sino(:,1,:));
dyn_range = [0 max( sino(:))];
imsc( im, dyn_range);
axis tight equal; colorbar; nn = nn + 1;
title( sprintf( 'sino(:,1,:)' ) )
name = sprintf( 'phase_map_n%u', 1);
save_png( im, name)
% retrieved phase map
set( figure(nn),  'Name',  'retrieved phase maps' );
im = squeeze(sino_retr(:,1,:)) ;
imsc( im );
axis tight equal; colorbar; nn = nn + 1;
title( sprintf( '[sino retr](:,1,:)' ) )
name = sprintf( 'phase_map_retrieved_n%u', 1);
save_png( im, name)

% intensity w/o noise
set( figure(nn),  'Name',  'intensity w/o noise' );
im = squeeze(int(:,1,:));
imsc( im );
axis tight equal; colorbar; nn = nn + 1;
title( sprintf( '[int](:,1,:)' ) )
name = sprintf( 'int_wo_noise_n%u', 1);
save_png( im, name)

% intensity with noise
set( figure(nn),  'Name',  'intensity with noise' );
im = squeeze(int_noise(:,1,:));
imsc( im );
axis tight equal; colorbar; nn = nn + 1;
title( sprintf( '[int noise prop](:,1,:)' ) )
name = sprintf( 'int_n%u', 1);
save_png( im, name)

% FT intensity
set( figure(nn),  'Name',  'FT of propagated projections' );
im = squeeze(int_noise(:,1,:));
im = log( 1 + abs( fftshift( fft2( SubtractMean( im ) ) ) ) );
imsc( im );
axis tight equal; colorbar; nn = nn + 1;
title( sprintf( 'FT[int noise](:,1,:)' ) )
name = sprintf( 'FT_propagated_projection_n%u', 1);
save_png( im, name)

% Filter
if ~strcmp( coherence, 'full' )
    set( figure(nn),  'Name',  'Fourier space blurring filter' );
    im = fftshift( Binning( filt ) / 4 );
    imsc( im );
    save_png( im, 'filter_partial_coherence');
    axis tight equal; colorbar; nn = nn + 1;
end

set( figure(nn),  'Name',  'reconstruction: difference maps |reco - phan|' );
sli = round( size( reco, 3) / 2 );
im = squeeze(diff_reco_phan(:,:,sli)) ./ max2( phan(:,:,sli) );
imsc( im );
axis tight equal; colorbar; nn = nn + 1;
title( sprintf( '1 / max(phan(:,:,%u)) * |reco - phan|(:,:,%u)', sli, sli ) )
name = sprintf( 'abs_tomo-phan_z%u_%2.0f%%ofMaxPhan', 1, 100*max( im(:)));
save_png( im, name)

% Sinograms
% % simulated phase sinogram
% set( figure(nn),  'Name',  'projection: sino' );
% sli = round( size( sino, 3) / 2 );
% imsc( squeeze(sino(:,:,sli)), [0 max( sino(:))]);
% axis tight equal; colorbar; nn = nn + 1;
% title( sprintf( 'sino(:,:,%u)', sli ) )
