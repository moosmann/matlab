function im = RotAxisSymmetricCropping(im, rot_axis_pos, dim)
% Crop image symmetrically around rotation axis without resampling such
% that the rotation axis is exactly centered within a pixel or between two
% pixels.
%
% ARGUMENTS
% im: array, 1D or 2D.
% rot_axis_pos: scalar. Position of the rotation axis.
% dim: 1 or 2. Default: 2. Dimension that should be cropped
%
% Written by Julian Moosmann, last version 2013-10-21. Modified: 2016-12-01

if nargin < 3
    dim = 2;
end

%% MAIN %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% round to next half integer
rot_axis_pos = round( 2* rot_axis_pos ) / 2;

num_pix = size( im, dim );
cen = num_pix / 2 + 0.5;
offset = rot_axis_pos - cen;

%% Rot axis is centered
if isequal( offset, 0 )    
    return

%% Rotation axis is off center
else
    % distance to closest edge
    d = min( rot_axis_pos - 1, num_pix - rot_axis_pos);

    %% Rot axis is right to center    
    if offset > 0
        start = rot_axis_pos - d;
        switch dim            
            case 1
                im = im(start:end, :);
            case 2
                im = im(:, start:end);
        end
    %% Rot axis is left to center        
    else
        end_new = rot_axis_pos + d;
        switch dim            
            case 1
                im = im(1:end_new, :);
            case 2
                im = im(:, 1:end_new);
        end        
    end
end
    
% INCONSISTENT AND THUS OBSOLET CODE
% % Rot axis is centered
% if rot_axis_pos == size( im, dim ) / 2
%     return
% 
% % Rot axis is left to center    
% elseif rot_axis_pos < size( im, dim ) / 2
%     if dim == 1
%         im = im(1:round(2*rot_axis_pos), :, :);
%     elseif dim == 2
%         im = im(:, 1:round(2*rot_axis_pos), :);
%     end
% % Rot axis is right to center
% else
%     newWidth = size( im, dim) - rot_axis_pos;
%     if dim == 1
%         im = im(end-round(2*newWidth)+1:end, :, :);        
%     elseif dim == 2
%         im = im( :, end-round(2*newWidth)+1:end, :);
%     end
% end