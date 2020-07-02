function f = iradonDesignFilter(filter_type, filter_length, freq_cutoff)
% Fourier Transform of the filter_type to be used to filter the projections
% before backprojection.
%
% INPUT ARGS:
% filter_type: string. Default: 'Ram-Lak'. Specifying the filter type.
% filter_length: scalar. Length of filter. Typically the detector width.
% freq_cutoff: scalar in [0,1]. Default: 1. Fraction of frequencies below
% the nyquist which we want to pass.
%
% OUTPUT ARGS: f: vector. Filter to use on the projections
%
% Taken from Matlab's iradon and modified by Julian Moosmann, 2016-09-29.
% Last modification: 2017-01-06
%
% f = iradonDesignFilter(filter_type, filter_length, freq_cutoff)

%% TODO: support odd number of pixels

%% Default arguments
if nargin < 1
    filter_type = 'Ram-Lak';
end
if nargin < 2
    filter_length = 1024;
end
if nargin < 3
    freq_cutoff = 1;
end

%% Main %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
order = 2 * floor( filter_length / 2);

filter_type = lower(filter_type);
if strcmpi(filter_type, 'none')
    f = ones(1, filter_length);
    return;
end

if ~strcmpi(filter_type, 'linear')
    % Create a bandlimited ramp filter (Eqn. 61 Chapter 3, Kak and Slaney)
    % - go up to the next highest power of 2.

    % 'order' is always even. 
    n = 0:(order / 2); 
    % 'filtImpResp' is the bandlimited ramp's impulse response (values for even n are 0)
    filtImpResp = zeros(1, (order / 2) + 1); 
    % Set the DC term 
    filtImpResp(1) = 1 / 4;
    % Set the values for odd n
    filtImpResp(2:2:end) = -1 ./ ((pi*n(2:2:end)).^2); 
    filtImpResp = [filtImpResp filtImpResp(end-1:-1:2)]; 

    f = 2 * real(fft(filtImpResp)); 
    f = f(1:(order/2)+1);
    %f = f / max(f(:));
else
    % Linear filter
    f = 0: order / 2;
    f = f / f(end);
end
% frequency axis up to Nyquist
w =  2 * pi * (0:size(f,2) - 1) / order;   

switch filter_type
    case {'ram-lak', 'linear'}
        % Do nothing    
    case 'shepp-logan'
        % be careful not to divide by 0:
        f(2:end) = f(2:end) .* (sin(w(2:end)/(2*freq_cutoff))./(w(2:end)/(2*freq_cutoff)));
    case 'cosine'
        f(2:end) = f(2:end) .* cos(w(2:end)/(2*freq_cutoff));
    case 'hamming'
        f(2:end) = f(2:end) .* (.54 + .46 * cos(w(2:end)/freq_cutoff));
    case 'hann'
        f(2:end) = f(2:end) .*(1+cos(w(2:end)./freq_cutoff)) / 2;
    otherwise
        error(message('images:iradon:invalidFilter'))
end

% Crop the frequency response
f( w > pi * freq_cutoff ) = 0;

% Symmetry of the filter
if mod( filter_length, 2 ) == 0
    f = [f' ; f(end-1:-1:2)'];
else
    f = [f' ; f(end:-1:2)'];
end

% Test case
%n=3;plot(-n:n-1,fftshift(iradonDesignFilter('Ram-Lak',2*n)),'.',-n:n,fftshift(iradonDesignFilter('Ram-Lak',2*n+1)),'-',-n:n,fftshift(iradonDesignFilter('linear',2*n+1)),'x')