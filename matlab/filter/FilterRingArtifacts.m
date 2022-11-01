function [imf, Vq, Vqf] = FilterRingArtifacts(im, rotAxisPos, rstride, thstride, verbose)

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
        rotAxisPos = [dim0/2, dim0/2];
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

figure( 'Name', 'polar image');
% subplot(1,2,1)
%title(sprintf('displacement offset: xq - x'))
imsc(rot90(Vq,1))
axis equal tight fill
xlabel('polar angle')
ylabel('radius')
xticks('auto'),yticks('auto')
drawnow

t = toc;

method = 'median';'guided';'wavelet-fft';
switch method
    case 'wavelet-fft'
        dec_levels = 3;
        wname = 'db25';
        sigma = 2.4;
        xc = round(size(Vq,1)/2);
        x1 = 1:xc;
        x2 = xc + x1;
        Vqf = cat(1, FilterStripesCombinedWaveletFFT( Vq(x1,:), dec_levels, wname, sigma ),FilterStripesCombinedWaveletFFT( Vq(x2,:), dec_levels, wname, sigma ));
    case 'guided'
        Vqf = imguidedfilter(Vq);
        
    case 'median'
        
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
        %
        % Inverse polar coordinates
        [thq, rq] = cart2pol(X, Y);
        
        % median filter width array
        med_filt_width = 2:10;
        for n = numel(med_filt_width):-1:1
            
            % median filter image
            fw = med_filt_width(n);
            Vm = medfilt2(Vq, [1, fw], 'symmetric');            
            polar_filt(:,:,n) = Vq - Vm;
            
            % Inverse polar transformation
            rings(:,:,n) = interp2(r, th, polar_filt(:,:,n), rq, thq, 'linear', 0);
            
            % Figure
            s = sprintf('filtered images: median kernel %u', fw);
            figure( 'Name', s, 'WindowState', 'maximized');
            subplot(1,2,1)
            title(sprintf('polar difference: width %u', fw))
            imsc(rot90(polar_filt(:,:,n),1))
            axis equal tight
            xlabel('polar angle')
            ylabel('radius')
            xticks('auto'),yticks('auto')
            
            subplot(1,2,2)
            title(sprintf('rings: width %u', fw))
            imsc(rings(:,:,n))
            axis equal tight
            xlabel('polar angle')
            ylabel('radius')
            xticks('auto'),yticks('auto')
        end
        
        drawnow
        
        %         delta = 0.01;
        %         mask = Vq./Vm > 1 + delta | Vq./Vm < 1 - delta;
        %
        %         Vq(mask) = Vm(mask);
        %
        %         Vqf = (Vq - 1) * (Vq_max - Vq_min) + Vq_min;
        
        % %% Radial median filtering
        % Vqf = medfilt2(Vq, [1, 15], 'symmetric');
end

% % Inverse polar transformation
% [thq, rq] = cart2pol(X, Y);
% imf = interp2(r, th, Vqf, rq, thq, 'linear', 0);

%% Figure filter image
figure( 'Name', 'filter image');

subplot(1,3,1)
title(sprintf('unfiltered image'))
imsc(im)
axis equal tight
xticks('auto'),yticks('auto')

subplot(1,3,2)
title(sprintf('filtered image'))
imsc(imf)
axis equal tight
xticks('auto'),yticks('auto')

subplot(1,3,3)
title(sprintf('rings'))
imsc(im - imf)
axis equal tight
xticks('auto'),yticks('auto')

drawnow

%% Print information to standard out
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
    domain(imf)
    
    fprintf('\n')
end

end
