function [imf, Vq, Vqf] = FilterRingsPolar(im, rot_axis_pos, rstride, thstride, verbose)

%% Default arguments %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Rotation axis position
if nargin < 2
    rot_axis_pos = [];
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
angle_offset = -pi/2;
%% MAIN %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% image dimensions
[dim0, dim1] = size(im);

% Rotation axis position
if isempty(rot_axis_pos)
    if dim0 ~= dim1
        rot_axis_pos = input('Image has non-square dimensions. Enter rotation axis vector:');
    else
        rot_axis_pos = [dim0/2, dim0/2];
    end
end
if isscalar(rot_axis_pos)
    if dim0 ~= dim1
        fprintf(' Rotation axis position is scalar but image is non-square.\n')
        returnq
    else
        rot_axis_pos = [rot_axis_pos, rot_axis_pos];
    end
end

% radial vector
rmax = sqrt(1) * max([dim0/2, dim1/2, rot_axis_pos, [dim0/2, dim1/2] - rot_axis_pos]);
rv = 0:rstride:rmax;
% angle vector
if thstride == -1
    thstride = 1 / rmax;
end
%thv = fliplr(pi:-thstride:-pi+thstride);
thv = fliplr(pi:-thstride:-pi);
% Original cartesian grid
[X, Y] = meshgrid((1:dim1) - rot_axis_pos(2), (1:dim0) - rot_axis_pos(1));
% Polar grid
[r, th] = meshgrid(rv, thv);
% Cartesian coordinates of the polar grid
[Xq, Yq] = pol2cart(th + angle_offset,r);
% Interpolate
tic;
Vq = interp2(X, Y, im, Xq, Yq, 'linear', 0);

figure( 'Name', 'polar image');
% subplot(1,2,1)
%title(sprintf('displacement offset: xq - x'))
imsc(rot90(FilterHisto(Vq),1))
axis equal tight fill
xlabel('polar angle')
ylabel('radius')
xticks('auto'),yticks('auto')
drawnow

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

%         Vq_min = min(Vq(:));
%         Vq_max = max(Vq(:));
%         Vq = 1 + normat(Vq);

% mask for two half rings
reg_par = 0.001;
a1 = floor(size(Vq,1)/2);
a2 = ceil(size(Vq,1)/2);
m1 = mean(Vq(1:a1,:),1);
m2 = mean(Vq(a2:end,:),1);

% Normalize to same mean
m1 = SubtractMean(m1) + 1;
m2 = SubtractMean(m2) + 1;

% Combined mask
m = cat(1,repmat(m1,a1,1),repmat(m2,a2,1));

% Inverse mask
mi = 1./(abs(m) + reg_par);
mi = (SubtractMean(mi) + 1)*mean2(Vq);

m1i = SubtractMean(1./(abs(m1) + reg_par)) + 1;
m2i = SubtractMean(1./(abs(m2) + reg_par)) + 1;
figure( 'Name', '1D mask and inverse');
plot([m1',m2',m1i',m2i'])
legend({'m1','m2','m1i','m2i'})


% Filter polar image
Vqf = Vq.*mi;
%Vqf = SubtractMean(Vqf) + mean(Vq);

% Inverse polar coordinates
[thq, rq] = cart2pol(X, Y);

% Inverse polar transformation
imf = interp2(r,th+angle_offset,Vqf,rq,thq+angle_offset,'linear',0);

z = imf == 0;
zm = mean2(imf(~z));
imf(z) = zm;

itool(FilterHisto(Vq'));itool(FilterHisto(mi'));
itool(FilterHisto(im));itool(FilterHisto(rot90(imf,1)));


%% Figure filter image
figure( 'Name', 'Image filtering');

subplot(1,3,1)
title(sprintf('unfiltered image'))
imsc(FilterHisto(im))
axis equal tight
xticks('auto'),yticks('auto')

subplot(1,3,2)
title(sprintf('filtered image'))
imsc(FilterHisto(imf))
axis equal tight
xticks('auto'),yticks('auto')

subplot(1,3,3)
title(sprintf('rings'))
imsc(FilterHisto(im - imf))
axis equal tight
xticks('auto'),yticks('auto')

drawnow

%% Print information to standard out
if verbose
    fprintf('Ring artifact filter:')
    fprintf('\n image size = [%g, %g]', dim0, dim1)
    fprintf('\n image center = [%g, %g]', dim0/2, dim1/2)
    fprintf('\n rotation axis position = [%g, %g]', rot_axis_pos)
    fprintf('\n radius vector: [min, max] = [%g, %g]', rv(1), rv(end))
    fprintf('\n number of angles = %g', length(thv))
    fprintf('\n angle stride = %g', thstride)
    fprintf('\n angle vector / degree: [min, max] = [%g, %g]', thv(1)*180/pi, thv(end)*180/pi)
    fprintf('\n size of query point matrices: Xq = [%g, %g], Yq = [%g, %g]', size(Xq), size(Yq))
    fprintf('\n elapsed time for interpolation: %g s', t)
    fprintf('\n Domain of original image: \n')
    domain(im)
    fprintf(' Domain of image in polar coordinates: \n')
    domain(imf)
    
    fprintf('\n')
end

end
