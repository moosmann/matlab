% Optimal distance.
dopt   = 0.5;
% Distance variation.
dd     = 0.03;
d1 = dopt-dd/2;
d2 = dopt+dd/2;
% Pixel size, energy
pixsiz = 0.7e-06;
energy = 10;
% Filter threshold.
bfthresh = 0.01;
% Create filter
f1 = SineFilter([2048 2048],[energy d1 pixsiz],bfthresh);
f2 = SineFilter([2048 2048],[energy dopt pixsiz],bfthresh);
f3 = SineFilter([2048 2048],[energy d2 pixsiz],bfthresh);
% Combine inverted filter. The more zeros the combination of inverted
% filters has the better the combination of these two is.
m12=((1-f1).*(1-f2));
m13=((1-f1).*(1-f3));
m23=((1-f2).*(1-f3));
% Compute mean of each filter and the combination.
fprintf('\nMean of inverted 1st filter at d=%4.3f: %.6f\n',d1,1-mean(f1(:)))
fprintf('Mean of inverted 2nd filter at d=%4.3f: %.6f\n',dopt,1-mean(f2(:)))
fprintf('Mean of inverted 3nd filter at d=%4.3f: %.6f\n',d2,1-mean(f3(:)))
fprintf('Mean of combination of inverted filter: %.6f (1st&3rd), %.6f (1st&2nd), %.6f (2n&3rd)\n',mean(m13(:)),mean(m12(:)),mean(m23(:)))
