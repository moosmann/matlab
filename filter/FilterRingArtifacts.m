function [out, Vq] = FilterRingArtifacts(im, rotAxisPos, rstride, thstride, verbose)

%% Default arguments %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Rotation axis position
if nargin < 2
    rotAxisPos = -1;
end
% Sampling width along radial direction
if nargin < 3
    rstride = 1;
end
% Sampling width along angular direction
if nargin < 4
    thstride = -1;
end
% Print information to standard out
if nargin < 5
    verbose = true;
end

%% MAIN %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% image dimensions
[dim0, dim1] = size(im);

% Rotation axis position
if rotAxisPos == -1
    if dim0 ~= dim1
        rotAxisPos = input('Image has non-square dimensions. Enter rotation axis vector:');
    else
        rotAxisPos = [dim0/2, dim0 / 2];
    end
end
if isscalar(rotAxisPos)
    if dim0 ~= dim1
        fprintf(' Rotation axis position is scalar but image is non-square.\n')
        return
    else
        rotAxisPos = [rotAxisPos, rotAxisPos];
    end
end

% radial vector
rmax = sqrt(1) * max([dim0/2, dim1/2, rotAxisPos, [dim0/2, dim1/2] - rotAxisPos]);
rv = 0:rstride:rmax;
% angle vector
if thstride == -1
    thstride = 1 / rmax;
end
thv = fliplr(pi:-thstride:-pi+thstride);
% Original cartesian grid
[X, Y] = meshgrid((1:dim1) - rotAxisPos(2), (1:dim0) - rotAxisPos(1));
% Polar grid
[r, th] = meshgrid(rv, thv);
% Cartesian coordinates of the polar grid
[Xq, Yq] = pol2cart(th, r);
% Interpolate
tic;
Vq = interp2(X, Y, im, Xq, Yq, 'linear', 0);
t = toc;

% Fourier transform polar image
% Vqf = fft2(Vq);
% Filtering
% Vqf(1,:) = squeeze(median([Vqf(1:2,:); Vqf(end,:)], 1));
% Vqf = real(ifft2(Vqf));

% Edge
% Vqf = Vq;
% edges = edge(Vq);
% Vqm = medfilt2(Vq, [3, 3], 'symmetric');
% Vqf(edges) = Vqm(edges);

Vq_min = min(Vq(:));
Vq_max = max(Vq(:));

Vq = 1 + normat(Vq);
Vm = medfilt2(Vq, [1, 6], 'symmetric');

delta = 0.01;
mask = Vq./Vm > 1 + delta | Vq./Vm < 1 - delta;

Vq(mask) = Vm(mask);

Vq = (Vq - 1) * (Vq_max - Vq_min) + Vq_min;

% %% Radial median filtering
% Vqf = medfilt2(Vq, [1, 15], 'symmetric');

% Inverse polar transformation
[thq, rq] = cart2pol(X, Y);
out = interp2(r, th, Vq, rq, thq, 'linear', 0);

% Print information to standard out
if verbose
    fprintf('Ring artifact filter:')
    fprintf('\n image size = [%g, %g]', dim0, dim1)
    fprintf('\n image center = [%g, %g]', dim0/2, dim1/2)
    fprintf('\n rotation axis position = [%g, %g]', rotAxisPos)
    fprintf('\n radius vector: [min, max] = [%g, %g]', rv(1), rv(end))
    fprintf('\n number of angles = %g', length(thv))
    fprintf('\n angle stride = %g', thstride)
    fprintf('\n angle vector / degree: [min, max] = [%g, %g]', thv(1)*180/pi, thv(end)*180/pi)
    fprintf('\n size of query point matrices: Xq = [%g, %g], Yq = [%g, %g]', size(Xq), size(Yq))
    fprintf('\n elapsed time for interpolation: %g s', t)
    fprintf('\n Domain of original image: \n')
    domain(im)
    fprintf(' Domain of image in polar coordinates: \n')
    domain(out)
    
    fprintf('\n')
end

end
