function [tlow, thigh, histo] = compression( vol, method, parameter, verbose )
% Return thresholds in order to compress, i.e. normalize, the dynamic range
% of volume array 'vol' according 'vol = (vol - tlow) ./ (thigh - tlow)'
% using 'method'. 
%
% ARGUMENTS
% vol : 3D-array
% method : string, default: 'full'. Method used to compress dynamic range:
%   'full' :  use the full dynamic range
%   'std' : 
% parameter : 1- or 2-component vector depending on method
% verbose : bool, default: 1. print information
%
% Written by J. Moosmann, 2017-06-14. Last version: 2017-06-15
%
% [tlow, thigh, histo] = compression( vol, method, parameter, verbose )

if nargin < 2
    method = 'outlier';
end
if nargin < 3
    parameter = [1 0];
end
if nargin < 4
    verbose = 1;
end

%% Main %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

histo = [];

% Compression method
switch lower( method )
    case 'full'
        % use the full dynamic range of the data
        tlow = min3( vol );
        thigh = max3( vol );
                
    case 'std'
        % use the dynamic range within several standard
        % deviations centered around the mean value of the data
        vol_mean = mean3( vol );
        vol_std = std3( vol );
        %num_std = parameter(1);
        tlow = vol_mean - num_stds * vol_std;
        thigh = vol_mean + num_stds * vol_std;        
    
    case 'threshold'
        % use dyanmic range within given lower and upper
        % thresholds
        tlow = parameter(1);
        thigh = parameter(2);
    case 'outlier'
        [~, tlow, thigh] = FilterOutlier( vol, [parameter(1), parameter(2)], '', 1, 0 );
    
    case 'histo'
        % crop dynamic range using histogram of volume ROI
        t2 = toc;
        PrintVerbose(verbose, '\nCompression by histogram thresholds:' )
        num_bins_max = 2048;
        bl = min( round( 1 / parameter(1) / 2 ) * 2, num_bins_max);
        bh = min( round( 1 / parameter(2) / 2 ) * 2, num_bins_max);
        num_bins = min( lcm( bl, bh ), num_bins_max);
        PrintVerbose(verbose, '\n lower/higher percentage : %g%%, %g%%', parameter*100 )
        PrintVerbose(verbose, '\n bins required for lower/higher percentage thresholds : %u, %u', bl, bh )
        PrintVerbose(verbose, '\n number of histogram bins (least common multiple) : %u', num_bins )
        [x, y, z] = size( vol );
        xx = round( x / 2 ) + (-ceil(0.9*x/2/sqrt(2)):floor(0.9*x/2/sqrt(2)));
        yy = round( y / 2 ) + (-ceil(0.9*y/2/sqrt(2)):floor(0.9*y/2/sqrt(2)));
        zz = round( z / 2 ) + (-ceil(z/4):floor(z/4));
        vol_roi_min = min3( vol(xx,yy,zz) );
        vol_roi_max = max3( vol(xx,yy,zz) );
        PrintVerbose(verbose, '\n volume ROI  min/max : %g, %g', vol_roi_min, vol_roi_max )
        [histo, edges] = histcounts( vol(xx,xx,zz), num_bins, 'BinLimits', [vol_roi_min vol_roi_max] );
        tlow = edges( round( num_bins / bl ) );
        thigh = edges( end - round( num_bins / bh ) );
        PrintVerbose(verbose, '\n done in %.0f s (%.2f min).',toc - t2, (toc - t2) / 60)
end

PrintVerbose(verbose, ' \n lower / higher threshold : %g, %g', tlow, thigh )

end
