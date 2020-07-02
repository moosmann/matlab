% Create images for presentation
%
% Requires reconstruction script 'p05_reco.m' to be exectued priorily as it
% assumes the existence of workspace variables

nn = 1;
data(nn).path = '/asap3/petra3/gpfs/p05/2016/commissioning/c20160803_001_pc_test/raw/phase_1000';
data(nn).slice = 6;657;
data(nn).proj.roi1 = 1:3056;
data(nn).proj.roi2 = 0*380+270:2250-0*231;
data(nn).flat.roi1 = 2301:2700;
data(nn).flat.roi2 = 821:1300;

nn = 1;
padding = 1;
bin = 2;

% Distances
dist_source_dcm = 50.9;
dist_source_sample = 82.7;
dist_dcm_sample = dist_source_sample - dist_source_dcm;

%% Real space
% Projections
slice = data(nn).slice;
roi1 = data(nn).proj.roi1;
roi2 = data(nn).proj.roi2;
froi1 = data(nn).flat.roi1;
froi2 = data(nn).flat.roi2;
p = rot90(proj(roi1,roi2,slice));
% Mean of best flat-file corrected projections
yshiftm = min( yshift, [], 2);
m = abs(yshiftm(:)) < 0.5;
fprintf( 'sum of proj-flat pairs with shift < 0.5: %u', sum(m))
pm = mean( proj, 3);
pm = rot90(pm(roi1,roi2));

% Flat fields
fm = mean( flat, 3);
fm = rot90( fm(roi1,roi2));
fmr = rot90( fm(froi2,froi1)); % !! fm already flipped

%% FFT
epsilon = 0;
fun = @(im) Binning( FilterHisto(fftshift(log( 10^(-epsilon) + abs( fft2( padarray( im , padding*size(im ), 'symmetric', 'post' ) ) ) ) ), 3), bin + 2*padding);
pf = fun( p );
pmf = fun( pm );
fmf = fun( fm );
fmrf = fun( fmr );

%% CTF binary mask
pixelsize = eff_pixel_size_binned*(1+padding);
[mask, mask_half] = BinaryMaskCTF( size(pf), energy, sample_detector_distance, pixelsize, regularization_parameter, binary_filter_threshold);

distance2 = sample_detector_distance + 0.75;
phase_shift = 0.2 * pi;
[mask2, mask_half2] = BinaryMaskCTF( size(pf), energy, distance2, pixelsize, regularization_parameter, binary_filter_threshold, phase_shift);

dist_dcm_detector = dist_dcm_sample + sample_detector_distance;
% Magnification factor
mag = (dist_dcm_detector + dist_source_dcm) / dist_source_dcm;
% Defocusing distance
defocus = dist_dcm_detector / mag;
distance3 = dist_dcm_detector;
pixelsize3 = pixelsize / mag;
[mask3, mask_half3] = BinaryMaskCTF( size(fmrf), energy, distance2, pixelsize, regularization_parameter, binary_filter_threshold, phase_shift);

%% Image plus half mask
pf_mh = normat(pf).*mask_half ;
pmf_mh = normat(pmf).*mask_half2;
fmf_mh = normat(fmf).*mask_half2;
fmrf_mh3 = normat(fmrf).*mask_half3;

%% Split images
pfm_f1 = [normat(pf(1:ceil(size(mask,1)/2),:)); normat(fmf(1+ceil(size(mask,1)/2):end,:))];
pfm_f2 = [normat(pf(:,1:ceil(size(mask,2)/2))), normat(fmf(:,1+ceil(size(mask,2)/2):end))];
pfm_f1r = pfm_f1.*mask;

%% Display
h1 = figure(1);

subplot(3,3,1)
imsc( p );
axis equal tight
title(sprintf('projection'))
colorbar

subplot(3,3,2)
imsc( pm );
axis equal tight
title(sprintf('mean projection'))
colorbar

subplot(3,3,3)
imsc( fm );
axis equal tight
title(sprintf('mean flat'))
colorbar

subplot(3,3,4)
imsc( pf );
axis equal tight
title(sprintf('FFT projection'))
colorbar

subplot(3,3,5)
imsc( pmf );
axis equal tight
title(sprintf('FFT mean projection'))
colorbar

subplot(3,3,6)
imsc( fmf );
axis equal tight
title(sprintf('FFT mean flat'))
colorbar

subplot(3,3,7)
imsc( pf_mh )
axis equal tight
title(sprintf('FFT projection * mask'))
colorbar

subplot(3,3,8)

% imsc(  pmf_mh )
% axis equal tight
% title(sprintf('FFT mean projection * mask'))
% colorbar
imsc(  fmrf_mh3 )
axis equal tight
title(sprintf('FFT mean flat roi * mask3'))
colorbar

subplot(3,3,9)
imsc( fmf_mh );
axis equal tight
title(sprintf('FFT mean flat * mask'))
colorbar

drawnow

%% Save images %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
im_path = sprintf('/home/moosmanj/images/%s/%s', beamtime_id, scan_name);
CheckAndMakePath( im_path )
%copyfile( logfile, im_path)

p_bin = Binning( p, bin);
write32bitTIF(sprintf('%s/proj.tif', im_path), p_bin);
imwritesc( p_bin, sprintf('%s/proj.png', im_path))

fm_bin = Binning( fm, bin);
write32bitTIF(sprintf('%s/flat_mean.tif', im_path), fm_bin);
imwritesc( fm_bin, sprintf('%s/flat_mean.png', im_path))

write32bitTIF(sprintf('%s/proj_fft.tif', im_path), pf);
imwritesc( pf, sprintf('%s/proj_fft.png', im_path))

im = fmf;
write32bitTIF(sprintf('%s/flat_mean_fft.tif', im_path), im);
imwritesc( im, sprintf('%s/flat_mean_fft.png', im_path))

im = pf_mh;
write32bitTIF(sprintf('%s/proj_fft_mask.tif', im_path), im);
imwritesc( im, sprintf('%s/proj_fft_mask.png', im_path))

im = fmf_mh;
write32bitTIF(sprintf('%s/flat_mean_fft_mask.tif', im_path), im);
imwritesc( im, sprintf('%s/flat_mean_fft_mask.png', im_path))

%% Parameter log %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
logfile_name = sprintf( '%s/parameter.log', im_path);
fid = fopen(logfile_name, 'w');
fprintf(fid, 'PROJECTIONS:\n');
fprintf(fid, 'energy : %s\n', energy);
fprintf(fid, 'pixel size : %s\n', eff_pixel_size);
fprintf(fid, 'distance : %s\n', sample_detector_distance);
fprintf(fid, '\nFLAT FIELDS\n');
fprintf(fid, 'distance : %s\n', distance2);
fprintf(fid, 'phase_shift : %s\n', phase_shift);
fprintf(fid, 'magnification : %s\n', mag);
fprintf(fid, 'defocus : %s\n', defocus);
fprintf(fid, '\nFOURIER TRANSFORMED IMAGES\n:');
fprintf(fid, 'function : %s\n', func2str(fun));
fprintf(fid, 'padding : %u\n', padding);
fprintf(fid, 'epsilon : %u\n', epsilon);
fprintf(fid, '\n');
fprintf(fid, 'scan_name : %s\n', scan_name);
fprintf(fid, 'beamtime_id : %s\n', beamtime_id);
fprintf(fid, 'scan_path : %s\n', scan_path);
fprintf(fid, 'camera : %s\n', cam);
fprintf(fid, 'num_dark_found : %u\n', num_dark);
fprintf(fid, 'num_ref_found : %u\n', num_ref_found);
fprintf(fid, 'num_ref_used : %u\n', num_ref_used);
fprintf(fid, 'num_proj_found : %u\n', num_proj_found);
fprintf(fid, 'num_proj_used : %u\n', num_proj_used);
fprintf(fid, 'raw_image_shape : %u %u\n', raw_im_shape);
fprintf(fid, 'raw_image_shape_binned : %u %u\n', raw_im_shape_binned);
fprintf(fid, 'binning_factor : %u\n', bin);
fprintf(fid, 'effective_pixel_size_mu : %g\n', eff_pixel_size * 1e6);
fprintf(fid, 'effective_pixel_size_binned_mu : %g\n', eff_pixel_size_binned * 1e6);
fprintf(fid, 'energy : %g eV\n', energy);
fprintf(fid, 'flat_field_correlation_area_1 : %u:%u:%u\n', flat_corr_area1(1), flat_corr_area1(2) - flat_corr_area1(1), flat_corr_area1(end));
fprintf(fid, 'flat_field_correlation_area_2 : %u:%u:%u\n', flat_corr_area2(1), flat_corr_area2(2) - flat_corr_area2(1), flat_corr_area2(end));
fprintf(fid, 'min_max_of_all_darks : %6g %6g\n', dark_min, dark_max);
fprintf(fid, 'min_max_of_median_dark : %6g %6g\n', dark_med_min, dark_med_max);
fprintf(fid, 'min_max_of_all_flats : %6g %6g\n', flat_min, flat_max);
fprintf(fid, 'min_max_of_all_corrected_flats : %6g %6g\n', flat_min2, flat_max2);
fprintf(fid, 'min_max_of_all_raws :  %6g %6g\n', raw_min, raw_max);
fprintf(fid, 'min_max_of_all_corrected_raws :  %6g %6g\n', raw_min2, raw_max2);
fprintf(fid, 'min_max_of_all_flat_corr_projs : %g %g \n', proj_min, proj_max);
% Phase retrieval
fprintf(fid, 'do_phase_retrieval : %u\n', do_phase_retrieval);
fprintf(fid, 'phase_retrieval_method : %s\n', phase_retrieval_method);
fprintf(fid, 'phase_retrieval_regularisation_parameter : %f\n', phase_retrieval_reg_par);
fprintf(fid, 'phase_retrieval_binary_filter_threshold : %f\n', phase_retrieval_bin_filt);
fprintf(fid, 'phase_padding : %u\n', phase_padding);
% Rotation
fprintf(fid, 'rotation_angle_full_rad : %f\n', rot_angle_full);
fprintf(fid, 'rotation_angle_offset_rad : %f\n', rot_angle_offset);
fprintf(fid, 'rotation_axis_offset_used : %f\n', rot_axis_offset);
fprintf(fid, 'rotation_axis_position_used : %f\n', rot_axis_pos);
fprintf(fid, 'raw_image_binned_center : %f\n', raw_im_shape_binned1 / 2);
fprintf(fid, 'rotation_correlation_area_1 : %u:%u:%u\n', rot_corr_area1(1), rot_corr_area1(2) - rot_corr_area1(1), rot_corr_area1(end));
fprintf(fid, 'rotation_correlation_area_2 : %u:%u:%u\n', rot_corr_area2(1), rot_corr_area2(2) - rot_corr_area2(1), rot_corr_area2(end));
fprintf(fid, 'date : %s', datetime);
fclose(fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf( '\n' )
