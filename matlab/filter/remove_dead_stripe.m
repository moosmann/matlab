function sinogram = remove_dead_stripe(sinogram, snr, sizeWin, residual, smooth_strength)
%REMOVE_DEAD_STRIPE  Remove unresponsive or fluctuating stripe artifacts in a sinogram.
%
% This is a MATLAB port of the Python function `remove_dead_stripe` from removal.py.
% Angular direction is assumed along rows (axis 0 in the Python code).
%
% Syntax:
%   S_out = remove_dead_stripe(S_in)
%   S_out = remove_dead_stripe(S_in, snr, sizeWin, residual, smooth_strength)
%
% Inputs (defaults match the Python code):
%   sinogram        - [nRow x nCol] single/double. Sinogram image.
%   snr             - scalar, default 3.0. Larger = less sensitive detection.
%   sizeWin         - integer, default 51. Median-window size (odd recommended).
%   residual        - logical, default true. If true, run residual cleanup
%                     with remove_large_stripe(...) at the end.
%   smooth_strength - integer, default 10. Window size of 1D uniform filter
%                     applied along rows (per column) to estimate local trend.
%
% Output:
%   sinogram        - Stripe-corrected sinogram.
%
% Notes:
%   - If you have a MATLAB function `detect_stripe(list_fact, snr)` on path,
%     it will be used. Otherwise, a robust fallback (median/MAD based) is used.
%   - The final residual cleanup calls remove_large_stripe(sinogram, snr, sizeWin).
%     If your MATLAB version of remove_large_stripe has a different signature,
%     adjust the call accordingly (see comment near the end).
%
% -------------------------------------------------------------------------

% ----- Defaults (mirror Python) -----
if nargin < 2 || isempty(snr),             snr             = 3.0;  end
if nargin < 3 || isempty(sizeWin),         sizeWin         = 51;   end
if nargin < 4 || isempty(residual),        residual        = true; end
if nargin < 5 || isempty(smooth_strength), smooth_strength = 10;   end

% Basic checks
validateattributes(sinogram, {'single','double'}, {'2d','nonempty'}, mfilename, 'sinogram', 1);
validateattributes(snr, {'single','double'}, {'scalar','>',1}, mfilename, 'snr', 2);
validateattributes(sizeWin, {'single','double'}, {'scalar','>=',1}, mfilename, 'sizeWin', 3);
validateattributes(smooth_strength, {'single','double'}, {'scalar','>=',1}, mfilename, 'smooth_strength', 5);

sinogram   = single(sinogram);                 % keep memory light (matches Python float32)
[nrow, ncol] = size(sinogram);

% ----- 1) Column-wise uniform smoothing along rows (like scipy.ndimage.uniform_filter1d) -----
% Reflect/symmetric boundary handling to mimic SciPy's default behavior.
if smooth_strength > 1
    k    = ones(smooth_strength, 1, 'single') / single(smooth_strength);
    pad  = floor(smooth_strength/2);
    sino_smooth = padarray(sinogram, [pad 0], 'symmetric', 'both');  % [nrow+2*pad x ncol]
    sino_smooth = conv2(sino_smooth, k, 'same');                    % back to [nrow x ncol]
    sino_smooth = sino_smooth(1+pad:end-pad,:);
else
    sino_smooth = sinogram;
end

% ----- 2) Column statistic & background via median filter -----
list_diff = sum(abs(sinogram - sino_smooth), 1);          % 1 x ncol
list_diff = single(list_diff);

% 1D median filtering (prefer movmedian; fallback to medfilt1 if needed)
if exist('movmedian','file') == 2
    list_diff_bck = movmedian(list_diff, max(1, round(sizeWin)));
elseif exist('medfilt1','file') == 2
    list_diff_bck = medfilt1(list_diff, max(1, round(sizeWin)));
else
    % Very simple fallback: median over a sliding window using conv (approximate)
    hw   = max(0, floor((round(sizeWin)-1)/2));
    list_diff_bck = list_diff;
    for c = 1:ncol
        lo = max(1, c-hw); hi = min(ncol, c+hw);
        list_diff_bck(c) = median(list_diff(lo:hi));
    end
end

nmean = mean(abs(list_diff_bck), 'all');
list_diff_bck(list_diff_bck == 0) = nmean;
list_fact = list_diff ./ list_diff_bck;                     % 1 x ncol

% ----- 3) Detect stripe columns (mask = 1 for "bad" columns) -----
if exist('detect_stripe','file') == 2
    list_mask = logical(detect_stripe(list_fact, snr));     % external if available
else
    list_mask = local_detect_stripe(list_fact, snr);        % robust fallback (below)
end

% Dilate by 1 neighbor on both sides (SciPy binary_dilation(iterations=1))
% Use convolution to avoid Image Processing Toolbox dependency
list_mask = conv(single(list_mask), single([1 1 1]), 'same') > 0;

% Clear first/last 2 columns (to ensure safe interpolation)
if ncol >= 2
    list_mask(1:min(2,ncol)) = false;
    list_mask(max(1,ncol-1):ncol) = false;
end

% ----- 4) Interpolate missing columns from the valid columns -----
x_valid = find(~list_mask);                     % good columns
x_miss  = find(list_mask);                      % bad columns
y_all   = (1:nrow).';                           % rows as column vector

% Only interpolate if there are some bad columns but fewer than ~1/3 of total
if ~isempty(x_miss) && numel(x_miss) < floor(ncol/3) && ~isempty(x_valid)
    % 2-D gridded interpolant on the subset of valid columns (linear in both dims)
    F = griddedInterpolant({y_all, x_valid}, sinogram(:, x_valid), 'linear', 'nearest');
    [YY, XX] = ndgrid(y_all, x_miss);
    sinogram(:, x_miss) = F(YY, XX);
end

% ----- 5) Optional residual cleanup (calls remove_large_stripe, as in Python) -----
if residual
    % Python: sinogram = remove_large_stripe(sinogram, snr, sizeWin)
    % If your MATLAB signature is remove_large_stripe(S, snr, size, drop_ratio, norm, ...),
    % use: sinogram = remove_large_stripe(sinogram, snr, sizeWin, 0.1, true);
    try
        sinogram = remove_large_stripe(sinogram, snr, sizeWin);
    catch
        % Fallback to a common MATLAB signature we used elsewhere:
        sinogram = remove_large_stripe(sinogram, snr, sizeWin, 0.1, true);
    end
end

end

% === Local robust fallback for stripe detection (median/MAD based) =========
function mask = local_detect_stripe(list_fact, snr)
% list_fact: 1 x n
% Returns logical mask: true where a stripe is detected.
    list_fact = double(list_fact(:).');          % row
    med  = median(list_fact);
    madv = median(abs(list_fact - med));
    sigma = 1.4826 * madv + eps;                 % robust scale (≈ std for Gaussian)
    mask = abs(list_fact - med) > (snr * sigma); % boolean row
end
