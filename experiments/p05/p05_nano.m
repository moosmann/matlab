do_phase = 1;
slice = 0;

int_path = '/asap3/petra3/gpfs/p05/2017/data/11003420/processed/nano1336_tomo_Ivory_A001_Zernike/tomopy_croped/normalized/radios';

if ~exist( 'p', 'var')
    p = permute( Readstack( int_path ) ,[ 2 1 3]);
end
proj = exp( p );

%% Phase retrieval filter
if do_phase
    padding = 1;
    method = 'tie';
    edp = [11000, 0.1 1e-6];
    reg_par = 0.5;
    bin_filt = 0.1;
    cutoff_frequ = 1;
    
    fprintf( '\n Phase retrieval:' );t = toc;    
    im_shape = [size(proj,1) , size(proj,2)];
    im_shape_pad = (1 + padding) * im_shape;
    [phase_filter, pha_appendix] = PhaseFilter( method, im_shape_pad, edp, reg_par, bin_filt, cutoff_frequ, 'single');
                
    parfor nn = 1:size(proj, 3)        
        im = padarray( proj(:,:,nn), padding * im_shape, 'symmetric', 'post' );
        im = -real( ifft2( phase_filter .* fft2( im ) ) );
        pha = im(1:im_shape(1), 1:im_shape(2));
        proj(:,:,nn) = pha;
    end
    figure( 'Name', 'phase')
    imsc(proj(:,:,1)')
    fprintf( ' done in %.2f min.', (toc - t) / 60)
    pause( 0.1)
end

%% Tomography

% Filter
fbp_filter_type = 'linear';
fbp_filter_padding = 1;
fpb_filter_freq_cutoff = 1;
fbp_filter_padding_method = 'symmetric';
take_neg_log = 0;

fprintf( '\n Filter sino:' );t = toc;
filt = iradonDesignFilter(fbp_filter_type, (1 + fbp_filter_padding) * size( proj, 1), fpb_filter_freq_cutoff);
proj_shape1 = size( proj, 1);
parfor nn =  1:size( proj, 2)
    im = proj(:,nn,:);
    im = padarray( NegLog(im, take_neg_log), fbp_filter_padding * [proj_shape1 0 0], fbp_filter_padding_method, 'post' );
    im = real( ifft( bsxfun(@times, fft( im, [], 1), filt), [], 1, 'symmetric') );
    proj(:,nn,:) = im(1:proj_shape1,:,:);
end
fprintf( ' done in %.2f min.', (toc - t) / 60)

% Backprojection
num_proj = size(proj, 3);
angles = pi * (0:num_proj - 1) / num_proj;
rot_axis_offset = -9.6124;
%[vol_shape, vol_size] = volshape_volsize( proj, [0.5 0.5 0.5], [-0.5 0.5 -0.5 0.5 -0.5 0.5], rot_axis_offset, 1);
[vol_shape, vol_size] = volshape_volsize( proj, [], [], rot_axis_offset, 1);

astra_pixel_size = 1;
link_data = 1;

fprintf( '\n Backproject:' );t = toc;
if slice > 0
    slice = round( size( proj, 2) * slice );
    sino = squeeze( proj(:,slice,:) );
    vshape = [vol_shape(1), vol_shape(2) 1];
    vsize = vol_size;
    vsize(5) = -0.5;
    vsize(6) = 0.5;
    vol = astra_parallel3D( sino, angles, rot_axis_offset, vshape,vsize, astra_pixel_size, link_data);
else
    vol = astra_parallel3D( permute(proj, [1 3 2]), angles, rot_axis_offset, vol_shape, vol_size, astra_pixel_size, link_data);
end
fprintf( ' done in %.2f min.', (toc - t) / 60)

%% Save
outpath = '/asap3/petra3/gpfs/p05/2017/data/11003420/processed/nano1336_tomo_Ivory_A001_Zernike/tomopy_croped/normalized/pseudo_phase_retrieval';
outpath = '/asap3/petra3/gpfs/p05/2017/data/11003420/processed/nano1336_tomo_Ivory_A001_Zernike/tomopy_croped/normalized/no_log';
CheckAndMakePath( outpath );
parfor nn= 1:size(vol,3)
    filename = sprintf('%s/reco_%06u.tif',outpath, nn);
    im = vol(:,:,nn);
    write32bitTIF( filename, im );
end

%% Display
fprintf( '\n Show:' );t = toc;
figure( 'Name', 'tomo')
imsc( squeeze( vol(:,:,ceil(size(vol,2)/2)) ) )
if size( vol, 3 ) > 1
    nimplayp( vol )
end
fprintf( '\nFinished')