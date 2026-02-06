function sinogram = remove_all_stripe(sinogram, snr, la_size, sm_size, drop_ratio, dim, varargin)
%REMOVE_ALL_STRIPE  Remove all types of stripe artifacts in a sinogram.
%
% This combines (in order):
%   1) remove_dead_stripe(...)   % algorithm 6
%   2) remove_large_stripe(...)  % algorithm 5
%   3) remove_stripe_based_sorting(...) % algorithm 3
%
% The angular direction is assumed along rows (axis 0 in the Python code).
% This function forwards optional filter "options" to the downstream calls.
%
% Parameters (with defaults matching the Python implementation):
%   sinogram   : 2D array (double/single). Sinogram image.
%   snr        : scalar, default 3.0.  Ratio (>1.0) for stripe detection;
%                larger is less sensitive.
%   la_size    : integer, default 51. Window size of the median filter used
%                to remove large stripes (recommend odd).
%   sm_size    : integer, default 21. Window size of the median filter used
%                to remove small-to-medium stripes (recommend odd).
%   drop_ratio : scalar in [0, 0.8], default 0.1. Fraction of rows to drop
%                when detecting large stripes to prevent false positives.
%   dim        : {1,2}, default 1. Dimension of the window for sorting-based
%                filtering (1 = columnwise, 2 = 2D window).
%
% Optional arguments (forwarded):
%   varargin   : name-value pairs or a struct of options to control the
%                smoothing filters in downstream methods. These are passed
%                directly to remove_large_stripe(...) and
%                remove_stripe_based_sorting(...).
%
% Returns:
%   sinogram   : 2D array. Stripe-removed sinogram.
%
% Notes:
%   - This is a direct MATLAB port of `remove_all_stripe` from removal.py
%     (N. T. Vo, 2021), preserving call order and defaults.  The three
%     helper functions must be present on your MATLAB path:
%       * remove_dead_stripe(sinogram, snr, size, residual, [smooth_strength])
%       * remove_large_stripe(sinogram, snr, size, drop_ratio, norm, varargin{:})
%       * remove_stripe_based_sorting(sinogram, size, dim, varargin{:})
%     Implementations of these are required for execution. [1](https://hereon-my.sharepoint.com/personal/julian_moosmann_hereon_de/Documents/Microsoft%20Copilot-Chatdateien/removal.py)
%
% Example:
%   S = remove_all_stripe(S, 3.0, 51, 21, 0.1, 1);
%
% -------------------------------------------------------------------------
% Default handling (to mirror Python defaults)
if nargin < 2 || isempty(snr),        snr        = 2.2;  end
if nargin < 3 || isempty(la_size),    la_size    = 251;   end
if nargin < 4 || isempty(sm_size),    sm_size    = 101;   end
if nargin < 5 || isempty(drop_ratio), drop_ratio = 0.1;  end
if nargin < 6 || isempty(dim),        dim        = 1;    end

% Basic input validation (lightweight)
validateattributes(sinogram, {'single','double'}, {'2d','nonempty'}, mfilename, 'sinogram', 1);
validateattributes(snr,        {'double','single'}, {'scalar','>',1},   mfilename, 'snr',        2);
validateattributes(la_size,    {'double','single'}, {'scalar','>=',1},  mfilename, 'la_size',    3);
validateattributes(sm_size,    {'double','single'}, {'scalar','>=',1},  mfilename, 'sm_size',    4);
validateattributes(drop_ratio, {'double','single'}, {'scalar','>=',0,'<=',0.8}, mfilename, 'drop_ratio', 5);
assert(ismember(dim, [1, 2]), 'remove_all_stripe:dim', 'dim must be 1 or 2.');

% -------------------------------------------------------------------------
% 1) Remove unresponsive/fluctuating stripes (algorithm 6)
%    Python: sinogram = remove_dead_stripe(sinogram, snr, la_size, residual=False)
sinogram = remove_dead_stripe(sinogram, snr, la_size, false);

% 2) Remove large stripes (algorithm 5)
%    Python: sinogram = remove_large_stripe(sinogram, snr, la_size, drop_ratio, **options)
%    In MATLAB we pass name-value/struct options via varargin.
sinogram = remove_large_stripe(sinogram, snr, la_size, drop_ratio, true, varargin{:});

% 3) Remove small-to-medium stripes via sorting (algorithm 3)
%    Python: sinogram = remove_stripe_based_sorting(sinogram, sm_size, dim, **options)
sinogram = remove_stripe_based_sorting(sinogram, sm_size, dim, varargin{:});

% Done
end
