function phan =  Ball(num_voxel,relative_origin,relative_radius)
% Create ball phantom

% Defaul arguments
if nargin < 1
    num_voxel = 100;
end
if nargin < 2
    relative_origin = 0.5;
end
if nargin < 3
    relative_radius = 0.2;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if numel( num_voxel ) == 1
    num_voxel = [num_voxel, num_voxel, num_voxel];
end
if numel( relative_origin ) == 1
    relative_origin = [ relative_origin, relative_origin, relative_origin ];
end

% Coordinates
[x, y, z] = meshgrid(1:num_voxel(2),1:num_voxel(1),1:num_voxel(3));

% Volume
phan = zeros(num_voxel);

% Radius
r = round( relative_radius * min( num_voxel(:) ) );

% Centre
x0 = round( relative_origin(1) * num_voxel(1) );
y0 = round( relative_origin(2) * num_voxel(2) );
z0 = round( relative_origin(3) * num_voxel(3) );

phan ( (x-x0).^2 + (y-y0).^2 + (z-z0).^2 - r^2 <= 0 ) = 1;