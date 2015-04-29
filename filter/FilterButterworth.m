function [g g1x g1y] = FilterButterworth(Dim_XY,CutoffFrequencies_XY,FilterOrder,ShowInfo)
% Butterworth filter

if nargin < 1
    Dim_XY = [100 200];
end
if nargin < 2 
    CutoffFrequencies_XY = [Dim_XY(1)/2 Dim_XY(2)/2];
end
if nargin < 3
    FilterOrder = 5;
end
if nargin < 4
    ShowInfo = 0;
end
% Image dimensions
dimx = Dim_XY(1);
dimy = Dim_XY(2);
% Cut-off frequencies
wcx = CutoffFrequencies_XY(1);
wcy = CutoffFrequencies_XY(2);
% Frequency range
w1x = 1:dimx;
w1y = 1:dimy;
% Filter gain
g1x = (1./(1+(w1x/wcx).^(2*FilterOrder)));
g1y = (1./(1+(w1y/wcy).^(2*FilterOrder)));
g1x = g1x.*fliplr(g1x);
g1y = g1y.*fliplr(g1y);
[gx gy] = meshgrid(g1y,g1x);
g = gx.*gy;
% Restore mean
g(1,1) = 1;
if ShowInfo
    % Plot filter.
    figure,plot(g1x);
    figure,plot(g1y);
    wx = [1 floor(dimx/2) wcx dimx];
    wy = [1 floor(dimy/2) wcy dimy];
    fprintf('Gain for wx = [%u %u %u %u]: [%f %f %f %f] or [%f %f %f %f] dB.\n',wx,gx(wx),10*log10(1./g1x(wx)))
    fprintf('Gain for wy = [%u %u %u %u]: [%f %f %f %f] or [%f %f %f %f] dB.\n',wy,gy(wy),10*log10(1./g1y(wy)))
    ishow(g)
end
