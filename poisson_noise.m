clc;
clear all;
% Read in standard MATLAB demo image.
originalImage = imread('cameraman.tif');
% Convert to double.
originalImage = double(originalImage);
subplot(1,3,1);
imshow(originalImage, []),colorbar;
max_value = max(max(originalImage)); % take the maximum value
normalizedImage = originalImage/ max_value; % normalize, now pixel values range from [0,1] (I know that there's no negative values in my image.
% The Poisson option of imnoise wants the values to be scaled by 1e-12,
% so let's have the values go from 0 to 10e-12:
normalizedImage = normalizedImage * 10e-12;
% Now we'd consider the image as a Poisson process with values of 0-10.
% Now plot.
subplot(1,3,2);
imshow(normalizedImage, []),colorbar;
noisyImage = 10e12* max_value * imnoise(normalizedImage, 'poisson'); % generate noisy image and scale back to the original intensities.
subplot(1,3,3);
imshow(noisyImage, []),colorbar;
