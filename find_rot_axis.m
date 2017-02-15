function [vol, sino] = find_rot_axis(proj, angles, offsets, slice, take_neg_log)

%% Default arguments
if nargin < 4
    slice = [];
end
if nargin < 5
    take_neg_log = 0;
end
num_pix = size(proj, 1);
subvol_shape = [num_pix, num_pix, 1];
subvol_size = [-num_pix/2 num_pix/2 -num_pix/2 num_pix/2 -0.5 0.5];
astra_pixel_size = 1;
link_data = 1;
rot_axis_tilt = 0;
number_of_stds = 3;
roi = 0.25;

%% Main
if isempty( slice )
    slice = round( size( proj, 2) / 2 );
end

sino = squeeze( proj(:, slice, :));

% Ramp filter
filt = iradonDesignFilter('Ram-Lak', num_pix, 0.9);

% Butterworth filter
[b, a] = butter(1, 0.5);
bw = freqz(b, a, numel(filt) );
filt = filt .* bw;

% Apply filters
sino = real( ifft( bsxfun(@times, fft( NegLog(sino, take_neg_log), [], 1), filt), [], 1, 'symmetric') );

vol = zeros(num_pix, num_pix, numel(offsets));
% Backprojection
for nn = 1:numel( offsets )
    offset = offsets(nn);
    vol(:,:,nn) = FilterHisto(astra_parallel3D( sino, angles, offset, subvol_shape, subvol_size, astra_pixel_size, link_data, rot_axis_tilt), number_of_stds, roi);
end
