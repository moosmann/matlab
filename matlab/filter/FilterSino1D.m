function [sino, G] = FilterSino1D(sino, alpha, verbose)

% Returns 1D filter for ring artifact suppression such that the filtered
% sinogram p is given by p = G * r, where r is the unfiltered sinogram and
% G a 1D filter, and * denotes convolution.

%% Default arguments
tic
if nargin < 2
    alpha = 0.002;
end
if nargin < 3
    verbose = true;
end

% Filter
vec_size = numel(sino) / 2;
G = alpha / sqrt(alpha * (4 + alpha)) * ( (2 + alpha - sqrt(alpha * (4 + alpha)) ) / 2 ).^abs(-vec_size:vec_size-1);
G = G(:);

% Convolution of Filter and sinogram
sino = fftshift(reshape(ifft( fft(G) .* fft(sino(:))), size(sino)));

% Print information
if verbose == true
    fprintf('\n elapsed time = %g s', toc)
    fprintf('\n size(G) = %g', size(G))
    fprintf('\n alpha = %g', alpha)
    fprintf('\n alpha([1, end]) = [%g, %g, %g]', G(1), G(vec_size+1), G(end))
    fprintf('\n')
end
