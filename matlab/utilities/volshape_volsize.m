function [vol_shape, vol_size] = volshape_volsize( proj, vol_shape, vol_size, rot_axis_offset, verbose, abs_rel_thresh)
% Returns shape and size of the volume to be reconstructed.
% 
% ARGUMENTS
% proj : projection array needed for image dimension
% vol_shape : empty or 3-element vector of relative or absolute numbers.
%   determines the number voxels of the reconstructed volume. 
% vol_size : empty or 6-element vector of relative of absolute numbers.
%   determines the extent of the reconstructed volume.
% rot_axis_offet : offset is taken into account to increase the volume to
%   be reconstructed, e.g. for scans with excentric rot axis without
%   stitching.
% abs_rel_thresh : threshold below which the numbers of vol_shape and
%   vol_size are interpreted as relative number.
% verbose : bool, default:0. print shape, size, and required memory
%
% Written by Julian Moosmann, 2018-01-10
%
% [vol_shape, vol_size] = volshape_volsize( proj, vol_shape, vol_size, rot_axis_offset, verbose, abs_rel_thresh)

%% Default arguments %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 4
    rot_axis_offset = [];
end
if nargin < 5
    verbose = 0;
end
if nargin < 6
    abs_rel_thresh = 10;
end

%% Main %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~isempty( vol_size ) && isempty( vol_shape )
    vol_shape = vol_size(2:2:end) - vol_size(1:2:end);
end

if isempty( rot_axis_offset )
    rot_axis_offset = 0;
end

%% Shape
vol_shape_hor = max( size( proj, 1), round( 2 * ( size( proj, 1)/2 + rot_axis_offset ) ) );
vol_shape_vert = size( proj, 2 );
if isempty( vol_shape )
    % default volume given by the detector width and height
    vol_shape = [vol_shape_hor, vol_shape_hor, vol_shape_vert];
else
    vol_shape_hor = size( proj, 1);
    if vol_shape(1) <=  abs_rel_thresh
        vol_shape(1) = round( 1 + vol_shape(1) * (vol_shape_hor - 1) );
    end
    if vol_shape(2) <=  abs_rel_thresh
        vol_shape(2) = round( 1 + vol_shape(2) * (vol_shape_hor - 1) );
    end
    if vol_shape(3) <=  abs_rel_thresh
        vol_shape(3) = round( 1 + vol_shape(3) * (vol_shape_vert - 1) );
    end
end

%% Size
if isempty( vol_size )    
%     vol_size(1) = ceil( -vol_shape(1)/2 );
%     vol_size(2) = floor( vol_shape(1)/2 );
%     vol_size(3) = ceil( -vol_shape(2)/2 );
%     vol_size(4) = floor( vol_shape(2)/2 );
%     vol_size(5) = ceil( -vol_shape(3)/2 );
%     vol_size(6) = floor( vol_shape(3)/2 );
    vol_size(1) = -vol_shape(1)/2;
    vol_size(2) = vol_shape(1)/2; 
    vol_size(3) = -vol_shape(2)/2;
    vol_size(4) = vol_shape(2)/2;
    vol_size(5) = -vol_shape(3)/2;
    vol_size(6) = vol_shape(3)/2;
else
    if abs( vol_size(1) ) < abs_rel_thresh
        vol_size(1) = sign( vol_size(1) ) * ( abs( vol_size(1) ) * (vol_shape_hor -1)  );
    end
    if abs( vol_size(2) ) < abs_rel_thresh
        vol_size(2) = sign( vol_size(2) ) * ( abs( vol_size(2) ) * (vol_shape_hor -1) );
    end
    if abs( vol_size(3) ) < abs_rel_thresh
        vol_size(3) = sign( vol_size(3) ) * ( abs( vol_size(3) ) * (vol_shape_hor -1)  );
    end
    if abs( vol_size(4) ) < abs_rel_thresh
        vol_size(4) = sign( vol_size(4) ) * ( abs( vol_size(4) ) * (vol_shape_hor -1) );
    end
    if abs(vol_size(5)) < abs_rel_thresh
        vol_size(5) = sign( vol_size(5) ) * max([( abs( vol_size(5) ) * (vol_shape_vert-1) ),0.5]);
    end
    if abs(vol_size(6))  < abs_rel_thresh
        vol_size(6) = sign( vol_size(6) ) * max([( abs( vol_size(6) ) * (vol_shape_vert-1) ),0.5]);
    end
end

%% Print info
PrintVerbose( verbose, '\n volume shape: [%g %g %g]', vol_shape )
PrintVerbose( verbose, '\n volume size: [%g %g %g %g %g %g]', vol_size )
PrintVerbose( verbose, '\n memory required before stitching and binning: %.2f GiB', prod( vol_shape ) * 4 / 1024^3 )
