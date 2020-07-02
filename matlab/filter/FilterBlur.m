function array = FilterBlur( array, median_neighborhood, disk_radius)
% Blur image.
%
% array : 1-, 2-, or 3D.
% median_neighboorhood : scalar or ndim-vector. values have to be odd for
%   3D arrays 
% disk_radius : scalar. Radius of simple disk filter.
% Written by Julian Moosmann, 2013-10-21. Mod: 2018-02-06

%% Defaults %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 2
    median_neighborhood = 3;
end
if nargin < 3
    disk_radius = 10;
end

%% MAIN %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Median
if sum( median_neighborhood ) > 0
    
    dim = ndims( array );
    if dim ~= numel( median_neighborhood )
        for nn = dim:-1:1
            neighborhood(nn) = median_neighborhood(1);
        end
    end
    
    
    switch dim
        case 1
            array = medfilt1(array, neighborhood, 'symmetric');
        case 2
            array = medfilt2(array, neighborhood, 'symmetric');
        case 3
            array = medfilt3(array, neighborhood, 'symmetric');
    end
end

% Disk
if disk_radius > 0
    array = imfilter( array, fspecial('disk',disk_radius), 'symmetric');
end
