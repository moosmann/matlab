function sinogram = remove_large_stripe(sinogram, snr, sizeWin, drop_ratio, normFlag, varargin)
%REMOVE_LARGE_STRIPE  Remove large stripe artifacts in a sinogram.
%
% MATLAB port of the Python function `remove_large_stripe` (N. T. Vo, 2021).
% The angular direction is assumed along rows (axis 0 in the Python code).
%
% Syntax:
%   S_out = remove_large_stripe(S_in)
%   S_out = remove_large_stripe(S_in, snr, sizeWin, drop_ratio, normFlag, ...)
%
% Inputs (defaults match the Python implementation):
%   sinogram   - [nRow x nCol] single/double. Sinogram image.
%   snr        - scalar, default 3.0. Larger = less sensitive stripe detection.
%   sizeWin    - integer, default 51. Median-filter window (odd recommended).
%   drop_ratio - scalar in [0, 0.8], default 0.1. Fraction of rows dropped
%                when estimating column factors to reduce false positives.
%   normFlag   - logical, default true. If true, normalize columns by the
%                detected factors before restoration.
%
% Optional smoothing options (name-value or struct):
%   - 'method','median_filter','para1',[1 size]       % default behavior
%   - 'method','gaussian_filter','para1',[1 size]     % Gaussian along rows
%   - 'method','uniform_filter','para1',[1 size]      % moving average along rows
%   - 'customFcn', @(A, varargin) ...                 % function handle to apply
%
% Returns:
%   sinogram   - Stripe-removed sinogram.
%
% Dependencies:
%   - Uses local helpers `sort_forward`, `sort_backward`, and a robust
%     fallback `local_detect_stripe` if your own `detect_stripe` is not on path.
%
% Reference:
%   This function reproduces the algorithmic steps of `remove_large_stripe`
%   in removal.py: sorting → smoothing → stripe detection → normalization →
%   inverse sorting and column replacement. [1](https://hereon-my.sharepoint.com/personal/julian_moosmann_hereon_de/Documents/Microsoft%20Copilot-Chatdateien/removal.py)
%
% -------------------------------------------------------------------------

% ----- Defaults (mirror Python) -----
if nargin < 2 || isempty(snr),        snr        = 3.0;  end
if nargin < 3 || isempty(sizeWin),    sizeWin    = 51;   end
if nargin < 4 || isempty(drop_ratio), drop_ratio = 0.1;  end
if nargin < 5 || isempty(normFlag),   normFlag   = true; end

% Basic checks
validateattributes(sinogram, {'single','double'}, {'2d','nonempty'}, mfilename, 'sinogram', 1);
validateattributes(snr,        {'double','single'}, {'scalar','>',1},   mfilename, 'snr',        2);
validateattributes(sizeWin,    {'double','single'}, {'scalar','>=',1},  mfilename, 'sizeWin',    3);
validateattributes(drop_ratio, {'double','single'}, {'scalar','>=',0,'<=',0.8}, mfilename, 'drop_ratio', 4);

sinogram = single(sinogram);                 % Python uses np.float32
[nrow, ncol] = size(sinogram);

% Clip drop_ratio and compute rows to drop (top & bottom)
drop_ratio = min(max(drop_ratio, 0.0), 0.8);
ndrop      = floor(0.5 * drop_ratio * nrow);

% ----- 1) Sort each column along rows (ascending) -----
[sino_sort, sino_index] = sort_forward(sinogram);   % like util.sort_forward(..., axis=0)

% ----- 2) Smooth the sorted sinogram (default: median along rows) -----
sino_smooth = sino_sort;
if isempty(varargin)
    % Median filter with a [sizeWin x 1] window along the row dimension.
    % movmedian along dim=1 achieves 1D median per column.
    if exist('movmedian','file') == 2
        w = max(1, round(sizeWin));
        sino_smooth = movmedian(sino_sort, w, 1);
    elseif exist('medfilt2','file') == 2
        sino_smooth = medfilt2(sino_sort, [sizeWin 1], 'symmetric');
    else
        % Minimal fallback: sliding median per column
        sino_smooth = slidingMedianCols(sino_sort, max(1, round(sizeWin)));
    end
else
    % Parse a simple options struct or name-value pairs
    opts = parseOptions(varargin{:});
    if isfield(opts, 'customFcn') && isa(opts.customFcn, 'function_handle')
        sino_smooth = opts.customFcn(sino_smooth, opts);
    elseif isfield(opts, 'method')
        method = lower(string(opts.method));
        para1  = []; if isfield(opts, 'para1'), para1 = opts.para1; end
        switch method
            case "median_filter"
                win = para1; if isempty(win), win = [1 sizeWin]; end
                if numel(win) == 2, win = win([2 1]); end  % expect [rows cols]
                if exist('medfilt2','file') == 2
                    sino_smooth = medfilt2(sino_smooth, win, 'symmetric');
                else
                    sino_smooth = movmedian(sino_smooth, max(1, win(1)), 1);
                end
            case "gaussian_filter"
                win = para1; if isempty(win), win = [1 sizeWin]; end
                r   = max(1, round(win(1)));                % rows kernel length
                sigma = max(1, round(win(end)))/6;          % heuristic sigma
                g = exp(-((-(r-1)/2:(r-1)/2).^2) / (2*sigma^2)); g = g(:)/sum(g);
                sino_smooth = conv2(sino_smooth, g, 'same');
            case "uniform_filter"
                win = para1; if isempty(win), win = [1 sizeWin]; end
                r   = max(1, round(win(1)));
                k   = ones(r,1,'single')/single(r);
                sino_smooth = conv2(sino_smooth, k, 'same');
            otherwise
                error('remove_large_stripe:options', ...
                      'Unknown method "%s". Use median_filter/gaussian_filter/uniform_filter or customFcn.', method);
        end
    else
        error('remove_large_stripe:options', 'Options must include "method" or "customFcn".');
    end
end

% ----- 3) Column factors from mid rows (exclude ndrop at top/bottom) -----
rowStart = ndrop + 1;
rowEnd   = nrow - ndrop;
if rowEnd < rowStart
    % If the data is too short, fall back to full-range means
    rowStart = 1; rowEnd = nrow;
end
list1 = mean(sino_sort(rowStart:rowEnd, :), 1);    % 1 x ncol
list2 = mean(sino_smooth(rowStart:rowEnd, :), 1);  % 1 x ncol

% Safe divide: ones where denominator is 0
list_fact = ones(1, ncol, 'like', list1);
nz        = (list2 ~= 0);
list_fact(nz) = list1(nz) ./ list2(nz);

% ----- 4) Detect stripe columns and dilate by one pixel -----
if exist('detect_stripe','file') == 2
    list_mask = logical(detect_stripe(list_fact, snr));     % external if available
else
    list_mask = local_detect_stripe(list_fact, snr);        % robust fallback (median/MAD)
end
list_mask = conv(single(list_mask), single([1 1 1]), 'same') > 0;  % binary_dilation(iterations=1)

% ----- 5) Optional normalization -----
if normFlag
    sinogram = sinogram ./ repmat(list_fact, nrow, 1);      % normalize columns
end

% ----- 6) Inverse sort & replace only the detected columns -----
sino_corr = sort_backward(sino_smooth, sino_index);          % like util.sort_backward(..., axis=0)
xlist_miss = find(list_mask > 0);
if ~isempty(xlist_miss)
    sinogram(:, xlist_miss) = sino_corr(:, xlist_miss);
end

end

% ===================== Helpers ===============================

function [s_sorted, idx] = sort_forward(S)
% Sort each column along rows; return sorted values and their original indices.
    [s_sorted, idx] = sort(S, 1, 'ascend');     % [nrow x ncol], [nrow x ncol]
end

function S = sort_backward(s_sorted, idx)
% Reconstruct original row order per column given the sorting indices.
    [nrow, ncol] = size(s_sorted);
    S = zeros(nrow, ncol, 'like', s_sorted);
    for j = 1:ncol
        S(idx(:,j), j) = s_sorted(:, j);
    end
end

function M = slidingMedianCols(Min, win)
% Sliding median along rows for each column (fallback without toolboxes).
    [nrow, ncol] = size(Min);
    M = zeros(nrow, ncol, 'like', Min);
    half = floor((win-1)/2);
    for j = 1:ncol
        col = Min(:, j);
        for i = 1:nrow
            lo = max(1, i-half); hi = min(nrow, i+half);
            M(i, j) = median(col(lo:hi));
        end
    end
end

function mask = local_detect_stripe(list_fact, snr)
% Robust stripe detector using median & MAD (≈ std for Gaussian).
    x = double(list_fact(:).');
    med  = median(x);
    madv = median(abs(x - med));
    sigma = 1.4826 * madv + eps;
    mask = abs(x - med) > (snr * sigma);
end

function opts = parseOptions(varargin)
% Parse options from name-value pairs or a single struct.
    if isscalar(varargin) && isstruct(varargin{1})
        opts = varargin{1};
        return;
    end
    if mod(numel(varargin), 2) ~= 0
        error('remove_large_stripe:options', 'Options must be name-value pairs or a struct.');
    end
    opts = struct();
    for k = 1:2:numel(varargin)
        name = string(varargin{k});
        val  = varargin{k+1};
        opts.(name) = val;
    end
end
