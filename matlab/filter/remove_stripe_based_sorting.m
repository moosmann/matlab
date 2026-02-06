function sinogram = remove_stripe_based_sorting(sinogram, sizeWin, dim, varargin)
%REMOVE_STRIPE_BASED_SORTING  Sorting-based stripe removal for sinograms.
%
% MATLAB port of the Python function `remove_stripe_based_sorting` (N. T. Vo).
% Steps:
%   1) Sort each column along rows (axis 0 in Python).
%   2) Smooth the sorted image (default: median; 1D along columns or full 2D).
%   3) Restore original row order via inverse sort.
%
% Syntax:
%   S_out = remove_stripe_based_sorting(S_in)
%   S_out = remove_stripe_based_sorting(S_in, sizeWin, dim, ...)
%
% Inputs (defaults match Python):
%   sinogram - [nRow x nCol] single/double. Sinogram image.
%   sizeWin  - integer, default 21. Median-filter window size (odd recommended).
%   dim      - {1,2}, default 1.
%              dim = 1 → 1D median along columns (window [1 x sizeWin]).
%              dim = 2 → 2D median (window [sizeWin x sizeWin]).
%
% Optional smoothing options (name-value pairs or struct), similar to Python **options:
%   'method','median_filter','para1',[r c]      % e.g., [1 sizeWin] (default behavior)
%   'method','gaussian_filter','para1',[sigma_row sigma_col]
%   'method','uniform_filter','para1',[r c]     % moving average
%   'customFcn',@(A,opts) ...                   % user-supplied function handle
%
% Returns:
%   sinogram - Stripe-smoothed sinogram after sort → filter → inverse sort.
%
% Notes:
%   - The function assumes stripes are vertical in the original sinogram.
%   - If Image Processing Toolbox functions are unavailable, robust fallbacks
%     are used (movmedian, separable Gaussian via conv1d, symmetric padding).
%
% -------------------------------------------------------------------------

% ----- Defaults -----
if nargin < 2 || isempty(sizeWin), sizeWin = 21; end
if nargin < 3 || isempty(dim),     dim     = 1;  end
assert(ismember(dim, [1 2]), 'dim must be 1 or 2.');

validateattributes(sinogram, {'single','double'}, {'2d','nonempty'}, mfilename, 'sinogram', 1);
validateattributes(sizeWin,  {'single','double'}, {'scalar','>=',1},  mfilename, 'sizeWin',  2);

% Work in single (Python uses np.float32)
S = single(sinogram);

% ----- 1) Sort each column along rows (ascending) -----
[Ssort, idx] = sort_forward(S);  % like util.sort_forward(..., axis=0)

% ----- 2) Smooth the sorted image -----
if isempty(varargin)
    % Default behavior: median filter
    if dim == 2
        % 2D median with window [sizeWin x sizeWin]
        if exist('medfilt2','file') == 2
            Ssort = medfilt2(Ssort, [sizeWin sizeWin], 'symmetric');
        else
            % Fallback: sequential 1D medians (approximate 2D median)
            Ssort = movmedian(Ssort, max(1, round(sizeWin)), 2);  % along columns
            Ssort = movmedian(Ssort, max(1, round(sizeWin)), 1);  % along rows
        end
    else
        % 1D median along columns (window [1 x sizeWin])
        if exist('medfilt2','file') == 2
            Ssort = medfilt2(Ssort, [1 sizeWin], 'symmetric');
        else
            Ssort = movmedian(Ssort, max(1, round(sizeWin)), 2);
        end
    end
else
    % Python-style options
    opts = parseOptions(varargin{:});
    if isfield(opts, 'customFcn') && isa(opts.customFcn, 'function_handle')
        Ssort = opts.customFcn(Ssort, opts);
    elseif isfield(opts, 'method')
        method = lower(string(opts.method));
        para1  = []; if isfield(opts, 'para1'), para1 = opts.para1; end
        switch method
            case "median_filter"
                % para1 = [rows cols]; default to [1 sizeWin] like Python
                win = para1; if isempty(win), win = [1 sizeWin]; end
                win = max(1, round(win));
                if exist('medfilt2','file') == 2
                    Ssort = medfilt2(Ssort, win([1 2]), 'symmetric');
                else
                    % Approximate: apply along dims indicated by window
                    if win(2) > 1, Ssort = movmedian(Ssort, win(2), 2); end % columns
                    if win(1) > 1, Ssort = movmedian(Ssort, win(1), 1); end % rows
                end

            case "gaussian_filter"
                % para1 = [sigma_row sigma_col]; default to [1 sizeWin]
                sig = para1; if isempty(sig), sig = [1 sizeWin]; end
                sig = double(sig(:).');
                if exist('imgaussfilt','file') == 2
                    % imgaussfilt allows [sigmaRow sigmaCol]
                    Ssort = imgaussfilt(Ssort, sig);
                else
                    % Separable Gaussian via conv1d with symmetric padding
                    if sig(1) > 0, Ssort = conv1d_symmetric(Ssort, gaussian1d(sig(1)), 1); end
                    if sig(2) > 0, Ssort = conv1d_symmetric(Ssort, gaussian1d(sig(2)), 2); end
                end

            case "uniform_filter"
                % para1 = [rows cols] kernel of ones / sum
                win = para1; if isempty(win), win = [1 sizeWin]; end
                win = max(1, round(win));
                Ssort = boxfilter_symmetric(Ssort, win);

            otherwise
                error('remove_stripe_based_sorting:options', ...
                      'Unknown method "%s". Use median_filter/gaussian_filter/uniform_filter or customFcn.', method);
        end
    else
        error('remove_stripe_based_sorting:options', 'Options must include "method" or "customFcn".');
    end
end

% ----- 3) Restore original row order via inverse sort -----
sinogram = sort_backward(Ssort, idx);

end

% ===================== Helpers ===============================

function [s_sorted, idx] = sort_forward(S)
% Sort each column along rows; return sorted values and original indices.
    [s_sorted, idx] = sort(S, 1, 'ascend');
end

function S = sort_backward(s_sorted, idx)
% Reconstruct original row order per column given sorting indices.
    [nrow, ncol] = size(s_sorted);
    S = zeros(nrow, ncol, 'like', s_sorted);
    for j = 1:ncol
        S(idx(:,j), j) = s_sorted(:, j);
    end
end

function opts = parseOptions(varargin)
% Parse options from name-value pairs or a single struct.
    if isscalar(varargin) && isstruct(varargin{1})
        opts = varargin{1}; return;
    end
    if mod(numel(varargin), 2) ~= 0
        error('remove_stripe_based_sorting:options', 'Options must be name-value pairs or a struct.');
    end
    opts = struct();
    for k = 1:2:numel(varargin)
        name = string(varargin{k});
        val  = varargin{k+1};
        opts.(name) = val;
    end
end

function g = gaussian1d(sigma)
% Normalized 1D Gaussian kernel (odd length ~ 6*sigma+1)
    L = max(3, 2*ceil(3*sigma)+1);
    x = -(L-1)/2:(L-1)/2;
    g = exp(-(x.^2) / (2*sigma^2));
    g = g(:) / sum(g);
end

function S = conv1d_symmetric(S, k, dim)
% 1D convolution along dimension dim (1=row, 2=col) with symmetric padding.
    half = floor(numel(k)/2);
    if dim == 1
        Spad = padarray(S, [half 0], 'symmetric', 'both');
        S    = conv2(Spad, k, 'valid');
    else
        Spad = padarray(S, [0 half], 'symmetric', 'both');
        S    = conv2(Spad, k', 'valid');
    end
end

function S = boxfilter_symmetric(S, win)
% Uniform (moving-average) filter with symmetric padding; win=[rows cols].
    r = max(1, round(win(1)));
    c = max(1, round(win(2)));
    kr = ones(r,1,'single')/single(r);
    kc = ones(1,c,'single')/single(c);
    S = conv1d_symmetric(S, kr, 1);
    S = conv1d_symmetric(S, kc', 2);
end
