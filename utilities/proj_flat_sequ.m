function [p, fig] = proj_flat_sequ( proj_roi, flat_roi, corr_matrix, proj_ind, num_flats, fig, dyn_range)
% Returns a 2D image which is the projection roi replicated and divided by
% the best matching flat-field rois.

if nargin < 4
    proj_ind = 1;
end
if nargin < 5
    num_flats = 10;
end
if nargin < 6
    fig = figure();
end
if nargin < 7
    dyn_range = [];
end

%% Main %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[corr_matrix_val, corr_matrix_pos] = sort( normat(corr_matrix), 2 );

pos = corr_matrix_pos(proj_ind,1:num_flats);
f = permute(flat_roi(:,:,pos), [2 1 3]);
f = f(:,:);
p = repmat( permute( proj_roi(:,:,proj_ind) ,[2 1 3] ), [1 num_flats] ) ./ f;


% Display
if nargin < 6
    fig = figure( 'Name', 'proj flat correlation');
end    
if isempty( dyn_range )
    dyn_range = [min(p(:)) max(p(:))];
end
imsc( p, dyn_range );
colorbar;
w = size( proj_roi, 1);
for nn = 1:num_flats
    t = text(5 + w * (nn -1 ), 25, sprintf( 'v %.0f', corr_matrix_val(proj_ind,nn) * 1000 ));
    t.Color = [1 1 1];
    t.FontSize = [16];
    t = text(5 + w * (nn -1 ), 50, sprintf( '# %u', corr_matrix_pos(proj_ind,nn) ));
    t.Color = [1 1 1];
    t.FontSize = [16];
end

t = text(5 , -25 + size( proj_roi,2) , sprintf( 'p # %u', proj_ind));
t.Color = [1 1 1];
t.FontSize = [16];


%truesize( fig )
