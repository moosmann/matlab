function [im, HotPixPercentage] = FilterHotPixel(im,FiltThreshHot_FiltThreshDark,printInfo,medianFilterRadius,filterDeadPixel)
%Deprecated. Calls FilterPixel. See FilterPixel for details.

[im, HotPixPercentage] = FilterPixel(im,FiltThreshHot_FiltThreshDark,printInfo,medianFilterRadius,filterDeadPixel);
