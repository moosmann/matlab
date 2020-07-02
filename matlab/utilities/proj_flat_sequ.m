% function [p, fig] = proj_flat_sequ( proj, flat, flat_corr_area1, flat_corr_area2, corr_matrix, proj_ind, num_flats, fig, dyn_range)
% % Returns a 2D image which is the projection roi replicated and divided by
% % the best matching flat-field rois.
%

proj_ind = 1200;
num_flats = 20;
corr_matrix = corr.mat;
    dyn_range = [0 1];

if ~exist( 'proj_ind', 'var' )
    proj_ind = 1;
end
if ~exist( 'num_flats', 'var' )
    num_flats = 6;
end
if ~exist( 'dyn_range', 'var' )
    dyn_range = [0 1];
end
if ~exist( 'corr_matrix', 'var' )
    corr_matrix = corr.mat;
end

%% Main %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[corr_matrix_val, corr_matrix_pos] = sort( normat(corr_matrix), 2 );

pos = corr_matrix_pos(proj_ind,1:num_flats);
f = permute(flat_roi(:,:,pos), [2 1 3]);

clear p
for nn = size( f, 3 ):-1:1
    p(:,:,nn) = permute( proj_roi(:,:,proj_ind) ,[2 1 3] ) ./ f(:,:,nn);
end

p = normarray( p );

%f = f(:,:);
%p = repmat( permute( proj_roi(:,:,proj_ind) ,[2 1 3] ), [1 num_flats] ) ./ f;

% Display
if ~exist( 'proj_flat_fig', 'var' )
    proj_flat_fig = figure( 'Name', 'proj flat correlation');
end
if isempty( dyn_range )
    dyn_range = [0.9*min(p(:)) 1.1*max(p(:))];
end
imsc( p(:,:), dyn_range );
axis tight
colorbar;
w = size( proj_roi, 1);
for nn = 1:num_flats
    t = text(5 + w * (nn -1 ), 25, sprintf( 'v %.0f', corr_matrix_val(proj_ind,nn) * 1000 ));
    t.Color = [1 1 1];
    t.FontSize = 16;
    t = text(5 + w * (nn -1 ), 50, sprintf( '# %u', corr_matrix_pos(proj_ind,nn) ));
    t.Color = [1 1 1];
    t.FontSize = 16;
end

% t = text(5 , -25 + size( proj_roi,2) , sprintf( 'p # %u, angle %g', proj_ind, angles(proj_ind)*180/pi));
% t.Color = [1 1 1];
% t.FontSize = 16;


%truesize( fig )
