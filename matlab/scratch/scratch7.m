filt = [ones( [size( proj, 1 ) 1] ); -ones( [size( proj, 1 ) 1] )];
filt(1) = 1;
%                filt = iradonDesignFilter(tomo.fbp_filter.type, (1 + tomo.fbp_filter.padding) * size( proj, 1), tomo.fbp_filter.freq_cutoff);
proj = dpc_phase;
proj_shape1 = size( proj, 1);
take_neg_log = 0;
padding = tomo.fbp_filter.padding;
padding_method = tomo.fbp_filter.padding_method;
parfor nn =  1:size( proj, 3)
    im = proj(:,:,nn);
    im = NegLog(im, take_neg_log);
    im = padarray( im, padding * [proj_shape1 0 0], padding_method, 'post' );
    im = fft( im, [], 1);
    im = bsxfun( @times, im, filt );
    im = real( ifft( im, [], 1, 'symmetric') );
    im = im(1:proj_shape1,:,:)
    proj(:,:,nn) = im;
end
fprintf( '\n finished\n' )
nimplayp( cat( 2, dpc_phase, proj ) )