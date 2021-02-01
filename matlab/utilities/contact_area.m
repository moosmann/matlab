function [c, m] = contact_area( vol, val1, val2 )
% Calculate the number of voxels (and surface) of the label with value
% 'val1' which are in contact with the labels with value 'val2'. Note that
% the number of contact voxels 'c' of 'contact_area( vol, val1, val2 )' is
% the same as 'contact_area( vol, val2, val1 )', the contact surface 'm' is
% not since the surface lies in label 1 in the former case and in the lable
% 2 in the latter case.
%
% ARGUMENTS
% vol : 3D array, segmented volume with distinct labels with values 'val1
%   and 'val2'
% val1 : scalar, value of label 1
% val2 : scalar, value of label 2
%
% OUTPUT
% c : scalar. Number of voxels of label 1 which are in contact with label 2 
% m : binary 3D array. Voxels of label 1 which are in contact with label 2 
%
% Written by Julian Moosmann
%
% function [c, m] = contact_area( vol, val1, val2 )

%% Default arguments %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 1
    %vol = zeros( [10 10 10] );
    %vol( 2:8, 2:8, 2:8 ) = 20;
    %vol( 4:6, 4:6, 4:6 ) = 30;
    vol = round( 100 * phantom3d( 'modified shepp-logan', 512 ) );
    val1 = 20;
    val2 = 30;
end

%% Main %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isequal( val1, val2 )
    error( 'use distinct values' )
end

tic
fprintf( 'Calculating contact surface:' )
m = zeros( size( vol ) );
[sx, sy, sz] = size( vol );
c = 0;
% Loop through all voxel
for x = 1:sx
    for y = 1:sy
        for z = 1:sz
            % voxel is in label 1
            if isequal( vol(x,y,z), val1 )
                if x+1 < sx && isequal( vol(x+1,y,z), val2 )
                    c = c + 1;
                    m(x,y,z) = m(x,y,z)+1;
                end
                if y+1 < sy && isequal( vol(x,y+1,z), val2 )
                    c = c + 1;
                   m(x,y,z) = m(x,y,z)+1;
                end
                if z+1 < sz && isequal( vol(x,y,z+1), val2 )
                    c = c + 1;
                    m(x,y,z) = m(x,y,z)+1;
                end
            end
            % voxel is in label 2
             if isequal( vol(x,y,z), val2 )
                if x+1 < sx && isequal( vol(x+1,y,z), val1 )
                    c = c + 1;
                    m(x+1,y,z) = m(x+1,y,z)+1;
                end
                if y+1 < sy && isequal( vol(x,y+1,z), val1 )
                    c = c + 1;
                   m(x,y+1,z) = m(x,y+1,z)+1;
                end
                if z+1 < sz && isequal( vol(x,y,z+1), val1 )
                    c = c + 1;
                    m(x,y,z+1) = m(x,y,z+1)+1;
                end
            end
        end
    end
end
m(m >= 1) = 1;
t = toc;
fprintf( '\n label 1 value : %g', val1 )
fprintf( '\n label 2 value : %g', val2 )
num_voxel = numel( vol );
fprintf( '\n number of voxel surfaces in contact : %u (%g%%)', c, 100 * c / num_voxel )
fprintf( '\n number of voxel: %u', num_voxel )
fprintf( '\n time elapsed : %f s', t )
fprintf( '\n' )
