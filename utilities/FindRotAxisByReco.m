function [rec, m] = FindRotAxisByReco( sino, angles, offset_range )
% Find rotation axis position via tomographic reconstruction
% 
% ARGUMENTS:
% sino: 2D sinogram, 3D sinograms not supported
% angles: 1D vector of angles corresponding to the sinogram
% offset_range: 1D vector, offset of the rotation axis position relative to the image center
% 
% RETURNS:
% rec: 3D array of reconstructed sinograms with different offsets
% m: 2D array, entropy and negativity metric for each offset
%
% Written Julian Moosmann, last version 2019-04-19
%
% [rec, m] = FindRotAxisByReco( sino, angles, offset_range )

%% Maina %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
astra_clear

% Volume geometry
vx = size( sino, 1 );
vy = vx;
vz = size( sino, 3 );

% Detector geometry
num_proj = numel( angles );
det_col_count = size( sino, 1 );
det_row_count = size( sino, 3 );
DetectorSpacingX = 1;
DetectorSpacingY = 1;

% Create ASTRA geometry vector
vectors = zeros( numel( angles ), 12);

rec = zeros( [vx vy numel( offset_range) ] , 'single' );
m = zeros( [2 numel( offset_range) ] , 'single' );

for oo = 1:numel( offset_range )
    
    
    offset = offset_range( oo );
    
    for nn = 1:num_proj
        
        theta = angles(nn);
        
        % ray direction
        vectors(nn,1) = sin( theta );
        vectors(nn,2) = -cos( theta );
        vectors(nn,3) = 0;
        
        % center of detector
        vectors(nn,4) = -offset * cos( theta );
        vectors(nn,5) = -offset * sin( theta );
        vectors(nn,6) = 0;
        
        % vector from detector pixel (0,0) to (0,1)
        vectors(nn,7) = cos( theta ) * DetectorSpacingX;
        vectors(nn,8) = sin( theta ) * DetectorSpacingX;
        vectors(nn,9) = 0;
        
        % vector from detector pixel (0,0) to (1,0)
        vectors(nn,10) = 0;
        vectors(nn,11) = 0;
        vectors(nn,12) = DetectorSpacingY;
        
    end
    
    % ASTRA projection geometry
    proj_geom = astra_create_proj_geom('parallel3d_vec',  det_row_count, det_col_count, vectors);
    
    % ASTRA volume geometry
    vol_geom = astra_create_vol_geom( vy, vx, vz );
    
    % ASTRA volume object
    vol_id = astra_mex_data3d('create', '-vol', vol_geom);
    
    sino_id = astra_create_sino3d_cuda(vol_id, proj_geom, vol_geom);
    
    % Filter sino
    pad_method = 'symmetric';'replicate';0;'none';
    sinof = FilterSinoForBackproj(sino, 1, 'Ram-Lak', pad_method, 'twice');
    astra_mex_data3d('set', sino_id, sinof)
    
    % ASTRA config struct
    cfg = astra_struct('BP3D_CUDA');
    cfg.ProjectionDataId = sino_id;
    cfg.ReconstructionDataId = vol_id;
    astra_mex_algorithm('create', cfg);
    
    % ASTRA create algorithm object from configuration struct
    bp_id = astra_mex_algorithm('create', cfg);
    
    % ASTRA backprojection
    astra_mex_algorithm('iterate', bp_id, 1);
    
    % Fetch data from ASTRA memory
    im = astra_mex_data3d('get_single', vol_id);
    rec(:,:,oo) = - im * pi/(2*length(theta));
        
    % ROI
    mask_rad = 0.95;
    mask_val = 0;
    roi = double( MaskingDisc( im, mask_rad, mask_val) ) * 2^16;
    
    % Metrics
    m(oo,1) = mean( roi( roi <= 0 ) );
    m(oo,2) = entropy( roi );
        
end

%% PLot
figure( 'Name', 'Metrics' )
Y(:,1) = normat( m(:,1) );
Y(:,2) = normat( m(:,2) );
plot( log( 1 + Y ) );
legend( 'negativity', 'entropy' )

%nimplay( rec )

%% Minima of metrics
[~,pn] = min( m(:,1) );
[~,pe] = min( m(:,2) );
fprintf( '\n negativity: pos = %u, offset = %f', pn, offset_range(pn) )
fprintf( '\n entropy   : pos = %u, offset = %f', pe, offset_range(pe) )
for nn = 1:numel( offset_range )
    fprintf( '\n %3u: %f', nn, offset_range(nn) )
end

fprintf( '\n' )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function sino = FilterSinoForBackproj(sino, direction, filt_type, pad_method, filt_len, filt_frequ_cutoff, butterworth)
% Filter to apply to sinogram before backprojection accounting for the
% change of coordinate system from cylindrical to cartesian. 
%
% sino: 2-or-3-dimensional array
% direction: scalar. Default: 1. Dimesnion along the filter should be
% applied
% filt_type: string. Default: 'Ram-Lak'. Options: 'none', 'Ram-Lak',
% 'Shepp-Logan', 'Cosine', 'Hamming', 'Hann'. See iradonDesingFilter for
% Details.
% pad_method: scalar, 'none', 'replicate', or 'symmetric'. Default: 0.
%
% Currently only even number of pixels in the filtering direction is
% possible. iradonDesignFilter does not yet support odd pixels.
%
% Written by Julian Moosmann
% First version: 2016-09-29. Last modification: 2016-09-30
%
% sino = FilterSinoForBackproj(sino, direction, filt_type)

%% TODO: improve documentantion
%% TODO: improve how padding parameters are provided

%% Default arguments %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 2
    direction = 1;
end
if nargin < 3
    filt_type = 'Ram-Lak';
end
if nargin < 4
    pad_method = 'symmetric';
end
if nargin < 5
    filt_len = 'twice';
end
if nargin < 6
    filt_frequ_cutoff = 1;
end
if nargin < 7
    butterworth = 0.5;
end

%% Main %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if direction == 0
    return
end

len = size( sino, direction); 

% Padding method
if strcmpi(pad_method, 'none')
    filt_len = len;
else
    % Padding size
    if strcmpi(filt_len, 'twice')
        filt_len = 2 * len;
    elseif strcmpi(filt_len, 'nextpow2')
        filt_len = max(64,2^nextpow2(2*len));
    else
        if filt_len < len
            error('Length of filter FILT_LEN must not be smaller than %u', len)
        end        
    end
    
    pad_vec = zeros([1, ndims(sino)] );
    pad_vec(direction) = filt_len - len;
    sino = padarray(sino, pad_vec,  pad_method, 'post');
end

% Design filter
filt = iradonDesignFilter(filt_type, filt_len, filt_frequ_cutoff);
filt_vec = ones([1, ndims(sino)] );
filt_vec(direction) = filt_len;
filt = reshape( filt, filt_vec);

% Butterworth filter
if butterworth > 0
    [b, a] = butter(6, butterworth);
    bw = freqz(b, a, numel(filt) );
    filt = filt .* bw;
end

% Fourier transform
sino = fft( sino, [], direction); % zero frequency is at first position 

% Filter multiplication
sino = bsxfun(@times, sino, filt);
    
% Inverse Fourier transform
sino = real( ifft( sino , [], direction, 'symmetric') );

% Crop to original size
if ~strcmpi(pad_method, 'none')
    switch direction
        case 1
            sino( len+1:end, :, :) = [];
        case 2
            sino( :, len+1:end, :) = [];
        case 3
            sino( :, :, len+1:end) = [];
    end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [im, m] = MaskingDisc( im, radial_fraction, value)
% Keep central ellipse within the 'radial_fraction' of the largest ellipse
% fitting within the rectangular image and set the region outside of this
% ellipse to the mean within the disc or to 'value'.
%
% im: 2D-matrix
% radial_fraction: scalar in [0 1]. default: 0.95. defines the radius of the
% disc to be keep as a fraction of the largest ellipsis that fits within the
% rectangular image
% value : scalar or []. default: []. if [] uses the mean of disc
%
% Written by Julian Moosmann, last version:2017-04-05
%
% im = MaskingDisc( im, radial_fraction, value) 

%% Default arguments %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 1
    im = rand( 10, 14 );
end
if nargin < 2
    radial_fraction = 0.95;
end
if nargin < 3
    value = [];
end

%% Main %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[x, y] = meshgrid(-size(im,2)/2:size(im,2)/2-1,-size(im,1)/2:size(im,1)/2-1);

x = x - (x(1,1)+x(1,end))/2;
y = y - (y(1,1)+y(end,1))/2;

m = sqrt((x/x(1,end)).^2 + (y/y(end,1)).^2);

m = m < radial_fraction;

if isempty(value)
    % set value to inner mean
    value = sum(im(:).*m(:))/sum(m(:));
end

im = m.*im + (1-m)*value;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function array = normat( array, min_max)
% Normalize matrix "array" to the dynamic range between min_max(1) and
% min_max(2), else between [0,1].
%
% Writtten by Julian Moosmann, last version: 2014-06-25

%% Define range
if nargin == 2
    min_val = min_max(1);
    max_val = min_max(2);
else
    min_val = min(array(:));
    max_val = max(array(:));
end

%% Normalize
array = ( array - min_val )/ ( max_val - min_val );
