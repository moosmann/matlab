function im = AddNoise(im,NoiseLevel,deltaNoise)

if nargin < 2
    NoiseLevel = 5000;
end
if nargin < 3
    deltaNoise = 0.01;
end

imMin = min(im(:));
imMax = max(im(:));
imMean = mean(im(:));

domain(im,1,'original  ')
im = 1 + deltaNoise* (im - imMean) / (imMax - imMin);
domain(im,1,'normalised')

im = 1/NoiseLevel*double(imnoise(uint16(NoiseLevel*im),'poisson'));

domain(im,1,'noised    ')

im = (im - 1) * (imMax - imMin) / deltaNoise + imMean;
%im (im<0) = 0;

domain(im,1,'noise rest');
