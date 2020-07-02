% function stat=hstparameter(name,n1,n2,nvue,X,Y,pixsize,corr,refon,offset,originx,originy,originz,sizex,sizey,sizez,offangle,varargin)
% create HST_SLAVE PARAMETER FILE
% the same file which would have been created by hst_master
% parameter file is n2.par created in n1
% name : name of output volume
% n1 : directory with the projections, resulting parameterfile and volume
% n2 : prefix of the projections
% nvue : number of vues 
% X : size of projections in horizontal direction
% Y : size of projections in vertical direction
% pixsize : size of the pixels in cm !!!
% corr : flatfield correction + logarithm to be done ('YES') or already done ('NO')
% refon : flatfield every refon images
% offset : position of rotation axis with respect to middle in pixels
% originx : first pixel in x (starting from 1)
% originy : first pixel in y (starting from 1)
% originz : first pixel in z (starting from 1)
% sizex : size of volume in x direction
% sizey : size of volume in y direction
% sizez : size of volume in z direction
% offangle : offset angle in degrees
% 
%       refname = varargin{1};
%	parfilename = varargin{2};
%	horaxis = varargin{3};
%       correct_spikes_threshold = varargin{4};
%
% origin: wolfgang, elodie, peter
 
function stat=hst_parameter(name,n1,n2,nvue,X,Y,pixsize,corr,refon,offset,originx,originy,originz,sizex,sizey,sizez,offangle,varargin)

global FT_CORRECT_SPIKES_THRESHOLD

switch nargin
case 17
	refname = 'refHST';
case 18
	refname = varargin{1};
case 19
	refname = varargin{1};
	parfilename = varargin{2};
case 20
	refname = varargin{1};
	parfilename = varargin{2};
	horaxis = varargin{3};
case 21
        refname = varargin{1};
	parfilename = varargin{2};
	horaxis = varargin{3};
        correct_spikes_threshold = varargin{4};
endswitch

% build parfilename
if not(exist('parfilename','var'))
	parfilename=sprintf('%s/%s.par',n1,n2);
elseif isempty(parfilename)
	parfilename=sprintf('%s/%s.par',n1,n2);
else
	parfilename=sprintf('%s/%s',n1,parfilename);
end

% build refname

if not(exist('refname','var'))
	refname='refHST';
elseif isempty(refname)
	refname='refHST';
end

if not(exist('horaxis','var'))
	horaxis = 0; 
elseif isempty(horaxis)
	horaxis = 0;
end

if isempty(FT_CORRECT_SPIKES_THRESHOLD)
	correct_spikes_threshold = [];
else
	correct_spikes_threshold = FT_CORRECT_SPIKES_THRESHOLD;
end

% correction spikes
if !exist('correct_spikes_threshold','var')
    correct_spikes_threshold = [];
end 
if isempty(correct_spikes_threshold)
    correct_spikes = 'NO';
    correct_spikes_threshold = 0.04;
else
    correct_spikes = 'YES';
end



fid=fopen(parfilename,'w');

if (fid==-1)
	disp('problems writing hst-parameterfile ! Check write permissions');
	return
end



%%%%ID 19 default values

first=0;
last=nvue-1;
number_length_varies= 'NO';
length_of_numerical_part=4;
file_postfix= '.edf';
flatfield_changing='YES';
num_image_1=X;
num_image_2=Y;
image_pixel_size_1=pixsize*10000;
image_pixel_size_2=pixsize*10000;

file_interval=1;
if horaxis
	rotation_vertical= 'NO';
else
	rotation_vertical= 'YES';
end

subtract_background= corr;
correct_flatfield= corr;
take_logarithm= corr;
not_available = 'N.A.';
if strcmp(corr,'YES')
	correct = 1;
else
	correct = 0;
end

background_file= [n1,'/dark.edf'];
flatfield_file= 'N.A.';
ff_prefix= [n1, '/', refname];
ff_num_first_image=0;
ff_num_last_image=nvue;
ff_file_interval=refon;

angle_between_projections = 180/nvue;
rotation_axis_position= X/2+offset;
output_sino = 'NO';
output_rec = 'YES';
start_voxel_1=originx;
start_voxel_2=originy;
start_voxel_3=originz;
end_voxel_1=originx+sizex-1;
end_voxel_2=originy+sizey-1;
end_voxel_3=originz+sizez-1;
oversampling_factor=4;
angle_offset=offangle;
% value ok both for alpha and beta
sino_megabytes=400;
cache_kilobytes=256; 
% provisional trick to avoid that hst on almond reconstructs 4 small slices simultaneously and 'forgets' to use all processors
% correction rings
correct_rings = 'NO';
correct_rings_nb = 8; % number of elements (|f|>0) that will not be filtered from sinograms 
% correction motion of axis (filename should contain drift expressed in pixels)
correct_axis = 'NO';
correct_axis_file = 'filename';
% padding at edges 'E' edge values or '0' 
padding = 'E';
% move axis to the center 'Y' or not 'N'
axis_to_the_center = 'N';
graph_out='NO';

if name(1)=='/'
	output_file = name;
else
	output_file=sprintf('%s/%s',n1,name);
end



% write out parameter file

fprintf(fid,'! HST_SLAVE PARAMETER FILE\n\n');
fprintf(fid,'! Parameters defining the projection file series\n');
fprintf(fid,'FILE_PREFIX = %s/%s\n',n1,n2);
fprintf(fid,'NUM_FIRST_IMAGE = %d ! No. of first projection file\n',first);
fprintf(fid,'NUM_LAST_IMAGE = %d ! No. of last projection file\n',last);
fprintf(fid,'NUMBER_LENGTH_VARIES = %s\n',number_length_varies);
fprintf(fid,'LENGTH_OF_NUMERICAL_PART = %d ! No. of characters\n',length_of_numerical_part);
fprintf(fid,'FILE_POSTFIX = %s\n',file_postfix);
fprintf(fid,'FILE_INTERVAL = %d ! Interval between input files\n',file_interval);

fprintf(fid,'\n! Parameters defining the projection file format\n');
fprintf(fid,'NUM_IMAGE_1 = %d ! Number of pixels horizontally\n',num_image_1);
fprintf(fid,'NUM_IMAGE_2 = %d ! Number of pixels vertically\n',num_image_2);
fprintf(fid,'IMAGE_PIXEL_SIZE_1 = %f ! Pixel size horizontally (microns)\n',image_pixel_size_1);
fprintf(fid,'IMAGE_PIXEL_SIZE_2 = %f ! Pixel size vertically\n',image_pixel_size_2);

fprintf(fid,'\n! Parameters defining background treatment\n');

if correct
	fprintf(fid,'SUBTRACT_BACKGROUND = %s ! Subtract background from data\n',subtract_background);
	fprintf(fid,'BACKGROUND_FILE = %s\n',background_file);
else
	fprintf(fid,'SUBTRACT_BACKGROUND = %s ! No background subtraction\n',subtract_background);
	fprintf(fid,'BACKGROUND_FILE = %s\n',not_available);
end
	                               
fprintf(fid,'\n! Parameters defining flat-field treatment\n');
if correct
	fprintf(fid,'CORRECT_FLATFIELD = %s ! Divide by flat-field image\n',correct_flatfield);
	fprintf(fid,'FLATFIELD_CHANGING = %s ! Series of flat-field files\n',flatfield_changing);
	fprintf(fid,'FLATFIELD_FILE = %s\n',flatfield_file);
	fprintf(fid,'FF_PREFIX = %s\n',ff_prefix);
	fprintf(fid,'FF_NUM_FIRST_IMAGE = %d ! No. of first flat-field file\n',ff_num_first_image);
	fprintf(fid,'FF_NUM_LAST_IMAGE = %d ! No. of last flat-field file\n',ff_num_last_image);
	fprintf(fid,'FF_NUMBER_LENGTH_VARIES = %s\n',number_length_varies);
	fprintf(fid,'FF_LENGTH_OF_NUMERICAL_PART = %d ! No. of characters\n',length_of_numerical_part);
	fprintf(fid,'FF_POSTFIX = %s\n',file_postfix);
	fprintf(fid,'FF_FILE_INTERVAL = %d ! Interval between flat-field files\n',ff_file_interval);
else
	fprintf(fid,'CORRECT_FLATFIELD = %s ! No flat-field correction\n',correct_flatfield);
	fprintf(fid,'FLATFIELD_CHANGING = %s\n',not_available);
	fprintf(fid,'FLATFIELD_FILE = %s\n',not_available);
	fprintf(fid,'FF_PREFIX = %s\n',not_available);
	fprintf(fid,'FF_NUM_FIRST_IMAGE = %s\n',not_available);
	fprintf(fid,'FF_NUM_LAST_IMAGE = %s\n',not_available);
	fprintf(fid,'FF_NUMBER_LENGTH_VARIES = %s\n',not_available);
	fprintf(fid,'FF_LENGTH_OF_NUMERICAL_PART = %s\n',not_available);
	fprintf(fid,'FF_POSTFIX = %s\n',not_available);
	fprintf(fid,'FF_FILE_INTERVAL = %s\n',not_available);
end

fprintf(fid,'\nTAKE_LOGARITHM = %s ! Take log of projection values\n',take_logarithm);

fprintf(fid,'\n! Parameters defining experiment\n');
fprintf(fid,'ANGLE_BETWEEN_PROJECTIONS = %f ! Increment angle in degrees\n',angle_between_projections);
fprintf(fid,'ROTATION_VERTICAL = %s\n',rotation_vertical);
fprintf(fid,'ROTATION_AXIS_POSITION = %f ! Position in pixels\n',rotation_axis_position);

fprintf(fid,'\n! Parameters defining reconstruction\n');
fprintf(fid,'OUTPUT_SINOGRAMS = %s ! Output sinograms to files or not\n',output_sino);
fprintf(fid,'OUTPUT_RECONSTRUCTION = %s ! Reconstruct and save or not\n',output_rec);
fprintf(fid,'START_VOXEL_1 =   %4d ! X-start of reconstruction volume\n',start_voxel_1);
fprintf(fid,'START_VOXEL_2 =   %4d ! Y-start of reconstruction volume\n',start_voxel_2);
fprintf(fid,'START_VOXEL_3 =   %4d ! Z-start of reconstruction volume\n',start_voxel_3);
fprintf(fid,'END_VOXEL_1 =   %4d ! X-end of reconstruction volume\n',end_voxel_1);
fprintf(fid,'END_VOXEL_2 =   %4d ! Y-end of reconstruction volume\n',end_voxel_2);
fprintf(fid,'END_VOXEL_3 =   %4d ! Z-end of reconstruction volume\n',end_voxel_3);
fprintf(fid,'OVERSAMPLING_FACTOR = %d ! 0 = Linear, 1 = Nearest pixel\n',oversampling_factor);
fprintf(fid,'ANGLE_OFFSET = %f ! Reconstruction rotation offset angle in degrees\n',angle_offset);
fprintf(fid,'CACHE_KILOBYTES = %d ! Size of processor cache (L2) per processor (KBytes)\n',cache_kilobytes);
fprintf(fid,'SINOGRAM_MEGABYTES = %d ! Maximum size of sinogram storage (megabytes)\n',sino_megabytes);

% adding extra features for PyHSTmpi
fprintf(fid,'\n! Parameters extra features PyHST\n');
fprintf(fid,'DO_CCD_FILTER = %s ! CCD filter (spikes)\n',correct_spikes);
fprintf(fid,'CCD_FILTER = "%s"\n','CCD_Filter');
fprintf(fid,'CCD_FILTER_PARA = {"threshold": %f }\n',correct_spikes_threshold);
fprintf(fid,'DO_SINO_FILTER = %s ! Sinogram filter (rings)\n',correct_rings);
fprintf(fid,'SINO_FILTER = "%s"\n','SINO_Filter');
fprintf(fid,'ar = Numeric.ones(%d,''f'')\n',num_image_1);
fprintf(fid,'ar[0]=0.0\n');
if (correct_rings_nb > 0)
	fprintf(fid,'ar[2:%d]=0.0\n',2*(correct_rings_nb+1));
end
fprintf(fid,'SINO_FILTER_PARA = {"FILTER": ar }\n');
fprintf(fid,'DO_AXIS_CORRECTION = %s ! Axis correction\n',correct_axis);
fprintf(fid,'AXIS_CORRECTION_FILE = ''%s''\n',correct_axis_file);
fprintf(fid,'OPTIONS= { ''padding'':''%s'' , ''axis_to_the_center'':''%s''} ! Padding and position axis\n',padding,axis_to_the_center);   

fprintf(fid,'\n! Parameters defining output file / format\n');
fprintf(fid,'OUTPUT_FILE = %s\n',output_file);

fprintf(fid,'\n! Reconstruction program options\n');
fprintf(fid,'DISPLAY_GRAPHICS = %s ! No images\n',graph_out);

fclose(fid);

unix(['chmod 666 ' parfilename]);  













