ca
name = 'syn027_Ti_12w_47R_z0900_t010ms_dpc';
outpath = '/asap3/petra3/gpfs/p07/2019/data/11006991/processed/images/';
filename = [outpath name '_roi.mat'];
if exist( 'rect', 'var' )
    fprintf( '\n crop' )
    imc = imcrop( f(:,:,1) , rect );
else   
    fprintf( '\n permute ' )
    f = permute( flat, [2 1 3] );
    fprintf( '\n normat' )
    f = normat( f );
    figure()
    fprintf( '\n crop' )
    [imc, rect] = imcrop( f(:,:,1) );
    rect = round( rect );
end
%%
x = rect(2) + (0:rect(4));
y = rect(1) + (0:rect(3));
fr = f(x,y,:);
fprintf( '\nsize imc: %u %u', size( imc ) )
fprintf( '\nsize f: %u %u %u', size( fr ) )
fprintf( '\n' )
domain( imc )
domain( fr(:,:,1) )
imsc( cat( 2, imc,  fr(:,:,1) ) )

%% Load cropped flat field
if ~exist( 'fr', 'var')
    load( filename );
end
%% Reshape
s1 = s(1);
s2 = s(2);
s3 = s(3);
f1 = fr(:,:,1:5:end);
f2 = fr(:,:,2:5:end);
f3 = fr(:,:,3:5:end);
f4 = fr(:,:,4:5:end);
f5 = fr(:,:,5:5:end);
nimplay( f1(:,:,2:end) - f1(:,:,1:end-1) )
%% %% Reshape
ca
s = size( fr );
f1 = reshape( fr, [s1, 5 * s2, s3 / 5 ]);
fprintf( '\n reshape : %u %u %u', size( f1 ) );
nimplay( f1 )
d = ( f1(:,:,2:end) - f1(:,:,1:end-1) );
nimplay( d )
nimplay(f5./f5(:,:,1))
%%
fprintf( '\n Writing')
for nn =1:20
    im = f5(:,:,nn) ./ f5(:,:,1);
    im = normat( im );
    filename = sprintf( '%s%s_flat_pos5Bypos5Im1_%06u.png', outpath, name, nn );
    imwrite( im, filename );
end
