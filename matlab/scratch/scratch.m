% p = reco_parameter;
% p.angles = 1:10;
% 
% testobject( p )
% disp( p.angles )
% function testobject( p )
%     disp( p.angles )
%     p.angles = p.angles + 10;
%     disp( p.angles )
% end
ca
clear all
raw_path = '/asap3/petra3/gpfs/p05/2020/data/11008993/raw/05_tomo_aorta240693_dpc_30keV';
%raw_path = '/asap3/petra3/gpfs/p05/2020/data/11008993/raw/23_tomo_aorta240693_dpc_Sqr10mu_30keV';
single_block = 0;
dpc_range = 30; % micron
eff_pix = 0.916;% l;micron
dpc_range_pix = round( dpc_range / eff_pix );

sr = dir( [raw_path filesep '*ref.tif'] );
if single_block
    sr1 = sr(1:15);
else
    sr1 = sr(1:15:end);
end
filename1 = [ raw_path filesep sr1(1).name ];
im1f = imread( filename1 );
filename2 =  [ raw_path filesep sr1(end).name ];
im2f = imread( filename2 );


if single_block
    fprintf( 'Crop images w.r.t. dpc shift' )
    im1f = im1f(:,1:end-dpc_range_pix);
    im2f = im2f(:,dpc_range_pix+1:end);
end


fprintf( '\n filename : %s', filename1 )
fprintf( '\n filename : %s', filename2 )
nimplay( cat( 3, im1f, im2f ), 0, [ 1 2 3], 'Full image' )


[optimizer, metric] = imregconfig('monomodal');

s = 0.45;
[x,y] = size( im1f );
x1 = round( s * x );
x2 = round( (1-s) * x);
y1 = round( s * y );
y2 = round( (1-s) * y );

im1 = im1f(x1:x2, y1:y2);
im2 = im2f(x1:x2, y1:y2);
nimplay( cat( 3, im1, im2 ), 0, [ 1 2 3], 'ROI' )

% itool( im1 )
% itool( im2 )
% itool( im1 - im2 )


% Translation: [1 0 0, 0 1 0, t_x t_y 0]
% Scale: [s_x 0 0, 0 s_y 0, 0 0 1]
% Shear: [1 sh_y 0, sh_x 1 0, 0 0 1]
% Rotation: [cos(t) sin(t) 0, -sin(t) cos(t) 0, 0 0 1]

%tform = imregtform( im1, im2, 'rigid', optimizer, metric);
tform = imregtform( im1, im2, 'translation', optimizer, metric);
angle_rad = asin( tform.T(1,2) ) / 2;
t1 = tform.T(3,1);
t2 = tform.T(3,2);
fprintf( '\n transfomratin matrix ; \n' )
disp( tform.T )

im1wfull = imwarp( im1, tform,  'OutputView', imref2d( size( im2) ) );
nimplay( cat( 3, im1wfull, im2 ), 1, [ 1 2 3], 'Registered full trafo ROI' )

im1wfullfull = imwarp( im1f, tform,  'OutputView', imref2d( size( im2f) ) );
nimplay( cat( 3, im1wfullfull, im2f ), 1, [ 1 2 3], 'Registered full trafo Full image' )

tform_norot = tform;
tform_norot.T(1,1) = 1;
tform_norot.T(2,2) = 1;
tform_norot.T(1,2) = 0;
tform_norot.T(2,1) = 0;
im1w_norot = imwarp( im1, tform_norot,  'OutputView', imref2d( size( im2) ) );
nimplay( cat( 3, im1w_norot, im2 ), 1, [ 1 2 3], 'Registered no rotation' )

fprintf( '\n angle : %g rad, %g degree', angle_rad, angle_rad * 180 / pi)
fprintf( '\n translation: %g, %g,', t1, t2 )

fprintf( '\n shift from rotation : %g, %g', [x, y] * sin( angle_rad ) ), 
fprintf( '\n' )
                            
                            
% % test volume
% if ~exist( 'v', 'var' )
%     v = uint16( floor( 2^16 *rand( [4000, 4000, 20]) ) );
%     v = padarray( v, [size(v,1) size(v,2)], 'symmetric', 'post' );
%     r = zeros( size( v ) );
% end
% 
% num_z = size( v, 3);
% gpu_arr = 1:gpuDeviceCount;
% 
% % number of available GPU
% num_gpu = numel( gpu_arr );
% % total amount of gpu task (using each GPU for multiple task)
% num_gpu_par = 2 * num_gpu;
% % number of for iteration: outer loop
% num_for = ceil( num_z / num_gpu_par );
% 
% % outer loop
% for mm = 1:num_for
%     fprintf( '\nfor loop interation: %u', mm )
%     
%     % par loop indices
%     nn_init = (mm-1)*num_gpu_par + 1;
%     nn_end = min( nn_init + num_gpu_par - 1, num_z );    
%     parfor_range = nn_init:nn_end;
%     
% fprintf( '\n' )
%     % inner parloop
%     for nn = parfor_range
%     %parfor nn = parfor_range
%         
%         gpu_index = mod( nn,  num_gpu ) + 1;
%         
%         %fprintf( '\nparfor loop index: %u', nn )
%         %fprintf( ', gpu index: %u', gpu_index )
%         %gpuDevice( gpu_index );
%         
% %         im = v(:,:,nn);
% %         im = FilterPixelGPU( im );
% %         r(:,:,nn) = im;
%         
% 
%         r(:,:,nn) = FilterPixelGPU( v(:,:,nn) );
% 
%         
%         
%         
%     end
%     
% end
% 
