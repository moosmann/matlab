function vol = PhaseRetrieval3D( vol, phase_method, EnergyDistancePixelsize, reg_par, bin_filt_thresh, padding, output_precision )
% Phase retrieval on reconstructed 3D volume of non-phase-retrieved
% intensity maps.

% Written by Julian Moosmann, 2017-05-24. Last modification:

%% Default arguments %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 2
    phase_method = 'tie';
end
if nargin < 3
    %EnergyDistancePixelsize = [30e3 1 1e-6];
    EnergyDistancePixelsize = [67e3 0.8 4*1.27e-6];
end
if nargin < 4
    reg_par = 2.5;
end
if nargin < 5
    bin_filt_thresh = 0.2;
end
if nargin < 6
    padding = 0;
end
if nargin < 7
    output_precision = class( vol );
end

%% Main %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic;
fprintf( '\n 3D phase filter:')

phase_filter = PhaseFilter3D( phase_method, (padding + 1) * size( vol ), EnergyDistancePixelsize, reg_par, bin_filt_thresh, output_precision);
fprintf( ' %.1f s', toc )

%vol = real( ifft( ifft( ifft( phase_filter .* fft( fft( fft( vol, [], 1), [], 2), [], 3), [], 3), [], 2), [], 1) );
if padding
    
    [x,y,z] = size( vol );
    
    t = toc;
        fprintf( '\n padding:');
    tmp = padarray(vol, padding * [x, y, z], 'symmetric', 'post');
    fprintf( ' %.2f min', (toc -t )/60 )
    
    t = toc;
    fprintf( '\n fftn:')
    tmp = fftn( tmp );
    fprintf( ' %.2f min', (toc -t )/60 )
    
    t = toc;
    fprintf( '\n phase filter multiplication:')
    tmp = phase_filter .* tmp;
    fprintf( ' %.2f min', (toc -t )/60 )
    
    t = toc;
    fprintf( '\n ifftn:')
    tmp = ifftn( tmp );
    fprintf( ' %.2f min', (toc -t )/60 )
    
    t = toc;
    fprintf( '\n real:')
    tmp = real( tmp );
    fprintf( ' %.2f min', (toc -t )/60 )
    
    t = toc;
    fprintf( '\n crop:')
    vol = tmp(1:x, 1:y, 1:z);
    fprintf( ' %.2f min', (toc -t )/60 )
else
    vol = real( ifftn( phase_filter .* fftn( vol ) ) );
end
fprintf('\n 3D phase retrieval done in %u s (%.2f min)\n', round(toc), (toc)/60 )
