## Copyright (C) 2007 P. Cloetens
##
## This program is free software; you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 2 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program; if not, write to the Free Software
## Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

## pyhst_parfile
##      stat=pyhst_parfile(direc,prefix,nvue,dim_1,dim_2,pixsize,corr,refon,offset,varargin)
##      creates parameter file for processing with PyHST
##      requires 9 fixed arguments and accepts options given as ['parameter',value]
##      the parameter name corresponds in most cases to the PyHST variable, but in lowercase
##
##      arguments:
##      argument 1: direc : path of the scan directory
##      argument 2: prefix : scan prefix
##      argument 3: nvue : number of angles used in the reconstruction
##      argument 4: num_image_1 : number of pixels in the projections along the first axis (horizontal)
##      argument 5: num_image_2 : number of pixels in the projections along the second axis (vertical)
##      argument 6: image_pixel_size : pixel size in microns
##      argument 7: corr : do flatfield correction + log or not: 'YES' or 'NO'
##      argument 8: refon : interval between the flatfield files: e.g. 100
##      argument 9: offset : position of rotation axis with respect to middle in pixels
##
##      options:
##      'parfilename'   : parameter file ( default: direc/prefix.par )
##      'start_voxel_1'       : first pixel in x (starting from 1) ( default 1 )
##      'start_voxel_2'       : first pixel in y (starting from 1) ( default 1 )
##      'start_voxel_3'       : first pixel in z (starting from 1) ( default 1 )
##      'end_voxel_1'         : last pixel in x (starting from 1) ( default num_image_1 )
##      'end_voxel_2'         : last pixel in y (starting from 1) ( default num_image_1 )
##      'end_voxel_3'         : last pixel in z (starting from 1) ( default num_image_2 )
##      'angle_offset'  : offset angle in degrees ( default 0 )
##      'output_file'   : name of the volume file ( default: direc/prefix.vol )
##      'ff_prefix'     : prefix of the flatfield images ( default: direc/refHST )
##      'background_file' : name of the dark file ( default: = direc/dark.edf )
##      'horaxis'       : rotation axis horizontal (=constant y) (1) 
##                                          or not (=constant x) (0) ( default 0 )
##      'correct_spikes_threshold' : threshold value for ccd filter 
##                                  when value is given, ccd filter will be applied with this parameter
##                                  ( default : parameter is set to 0.04, but DO_CCD_FILTER = 'NO' )
##      'ccd_filter'    : filter for correction of spikes ( default: 'CCD_Filter' )
##      'ccd_filter_para' : parameters for ccd filter
##          when given overrides choice with correct_spikes_threshold
##      'do_axis_correction' : correct for sample motion ( default: 'NO' )
##      'axis_correction_file' : name of file with sample movement ( default: 'correct.txt' )
##                              first column: x motion [ , second column: y motion ]
##      'do_sino_filter': correction of rings by filtering the sinogram ( default: 'NO' )
##      'sino_filter'   : filter for correction of rings ( default: 'SINO_Filter' )
##      'correct_rings_nb' : number of elements (|f|>0) that will not be filtered from sinograms
##                          ( default 8 )
##      'padding'       : padding at edges with edge values 'E' or zero '0' ( default: 'E' )
##      'axis_to_the_center' : move axis to the center 'Y' or not 'N' (default: 'Y') !!!
##      'angle_between_projections' : angle between successive projections ( default 180/nvue )
##      'doubleffcorrection' : name of the file (e.g. filtered_mean.edf) for double flatfield
##                              subtracted after normal flatfield and logarithm
##      'do_projection_median' : calculate median of projections (after flat/log) ( default: 'NO' )
##      'projection_median_filename' : name of output edf-file with median ( default: median.edf )
##      'do_projection_mean' : calculate mean of projections (after flat/log) ( default: 'NO' )
##      'projection_mean_filename' : name of output edf-file with mean ( default: mean.edf )
##                              PyHST uses for both median and mean PROJECTION_MEDIAN_FILENAME
##
##      less used options: (= do not change, except if you know what you are doing)
##      'zeroclipvalue' : minimum value of radiographs after flatfield / before logarithm
##                          values below are clipped ( default 1e-9 )
##      'oneclipvalue' : maximum value of radiographs after flatfield / before logarithm
##                          values above are clipped ( default: not applied )
##      'num_first_image' : index of first image used in reconstruction ( default 0)
##      'num_last_image' : index of last image ( default nvue-1 )
##      'number_length_varies' : width index varies or not ( default: 'NO' )
##      'length_of_numerical_part' : width of image index ( default 4 )
##      'file_postfix'  : extension of the image files ( default: '.edf')
##      'file_interval' : downsample in angle ( default 1 )
##      'subtract_background' : subtract dark file 'YES' or not 'NO' ( default: corr )
##      'correct_flatfield' : divide by flatfield image 'YES' or not 'NO' ( default: corr )
##      'take_logarithm'  : use -log 'Y' or 'N' ( default: corr )
##      'output_sinograms' : output sinogram files or not ( default: 'NO' )
##      'output_reconstruction' : do the reconstruction or not ( default: 'YES' )
##      'oversampling_factor' : oversampling of projections before backprojection
##                              (0 = linear interpolation, 1 = nearest pixel) ( default 4 )
##      'sinogram_megabytes' : maximum memory in MB used for reading the sinograms ( default 400 )
##      'cache_kilobytes'   : obsolete ( default 256 ) 
##      'display_graphics'  : NOT IMPLEMENTED IN PYHST ( default: 'NO' )
##
##      See also : hst_parameter (old version), parcut
      
## Author: P. Cloetens <cloetens@esrf.fr>
##
## 2007-03-04 P. Cloetens <cloetens@esrf.fr>
## * Initial revision
## 2007-03-09 PC
## * use start_voxel instead of origin
## * use end_voxel instead of size
## 2007-03-12 PC
## * change start comment from ! (Fortran) to # (Python and Octave)
## * fopen as binary to avoid warning
## * add options ZEROCLIPVALUE / ONECLIPVALUE
## 2007-05-08 PC
## * do not read correct_spikes_threshold from FT_CORRECT_SPIKES_THRESHOLD
## * should be done in function that calls pyhst_parfile
## * correct errors in setting of do_ccd_filter
## * foresee case correct_spikes_threshold = Inf
## 2007-05-11 PC
## * correct bugs in case horaxis found by Paul Tafforeau and Christian Nemoz
 
function stat=pyhst_parfile(direc,prefix,nvue,num_image_1,num_image_2,image_pixel_size,corr,refon,offset,varargin)

correct_spikes_threshold_default = 0.04;

stat = 1;
fixedvariables = 9;
if nargin < fixedvariables
    help pyhst_parfile
    return
endif

###############################
##### ESRF default values #####
###############################
not_available = 'N.A.';

parfilename=sprintf('%s/%s.par',direc,prefix);
num_first_image = 0;
num_last_image = nvue-1;
number_length_varies = 'NO';
length_of_numerical_part = 4;
file_postfix = '.edf';
image_pixel_size_1=image_pixel_size_2=image_pixel_size;
file_interval=1;
subtract_background = correct_flatfield = take_logarithm = corr;
correct = strcmp(corr,'YES');
background_file = [direc,'/dark.edf'];
flatfield_changing ='YES';
flatfield_file = not_available;
ff_prefix= [direc, '/refHST'];
ff_num_first_image=0;
ff_num_last_image=nvue;
ff_file_interval=refon;
doubleffcorrection = [];
zeroclipvalue = [];
oneclipvalue = [];
angle_between_projections = 180/nvue;
horaxis = 0;
output_sinograms = 'NO';
output_reconstruction = 'YES';
start_voxel_1 = start_voxel_2 = start_voxel_3 = 1;
end_voxel_1 = end_voxel_2 = num_image_1;
end_voxel_3 = num_image_2;
angle_offset = 0;
oversampling_factor = 4;
# limit the memory taken by each reconstruction process to ~ sinogram_megabytes
sinogram_megabytes=400;
# provisional trick to avoid that hst on almond reconstructs 4 small slices simultaneously and 'forgets' to use all processors
cache_kilobytes=256; 
# correction spikes
ccd_filter = 'CCD_Filter';
# correction rings
do_sino_filter = 'NO';
correct_rings_nb = 8; % number of elements (|f|>0) that will not be filtered from sinograms
sino_filter = 'SINO_Filter';
# correction motion of axis (filename should contain drift expressed in pixels)
do_axis_correction = 'NO';
axis_correction_file = 'correct.txt';
# padding at edges 'E' edge values or '0' 
padding = 'E';
# move axis to the center 'Y' or not 'N'
axis_to_the_center = 'Y';
do_projection_median = 'NO';
projection_median_filename = 'median.edf';
do_projection_mean = 'NO';
projection_mean_filename = 'mean.edf';
display_graphics='NO';
output_file = [direc '/' prefix '.vol'];

#######################################################
##### possibility to overwrite the default values #####
#######################################################

stat = parse_variables(varargin{:});


###################################################
##### check and calculate dependent variables #####
###################################################

if horaxis
    rotation_vertical= 'NO';
    rotation_axis_position= num_image_2/2+offset;
    end_voxel_12_limit = num_image_2;
    end_voxel_3_limit = num_image_1;
else
    rotation_vertical= 'YES';
    rotation_axis_position= num_image_1/2+offset;
    end_voxel_12_limit = num_image_1;
    end_voxel_3_limit = num_image_2;
end

# start_voxel should be > 0 ; end_voxel should be <= image dimension
if (start_voxel_1 < 1)
    printf('START_VOXEL_1 can not be zero or negative, forcing to 1\n')
    start_voxel_1 = 1;
endif
if (start_voxel_2 < 1)
    printf('START_VOXEL_2 can not be zero or negative, forcing to 1\n')
    start_voxel_2 = 1;
endif
if (start_voxel_3 < 1)
    printf('START_VOXEL_3 can not be zero or negative, forcing to 1\n')
    start_voxel_3 = 1;
endif
if (end_voxel_1 > end_voxel_12_limit)
    printf('END_VOXEL_1 can not be larger than image dimension, forcing to maximum\n')
    end_voxel_1 = end_voxel_12_limit;
endif
if (end_voxel_2 > end_voxel_12_limit)
    printf('END_VOXEL_2 can not be larger than image dimension, forcing to maximum\n')
    end_voxel_2 = end_voxel_12_limit;
endif
if (end_voxel_3 > end_voxel_3_limit)
    printf('END_VOXEL_3 can not be larger than image dimension, forcing to maximum\n')
    end_voxel_3 = end_voxel_3_limit;
endif

if (output_file(1)!='/')
    output_file=sprintf('%s/%s',direc,output_file);
endif
if (parfilename(1)!='/')
    parfilename=sprintf('%s/%s',direc,parfilename);
endif
if (background_file(1)!='/')
    background_file=sprintf('%s/%s',direc,background_file);
endif
if (ff_prefix(1)!='/')
    ff_prefix=sprintf('%s/%s',direc,ff_prefix);
endif

#############################
##### correction spikes #####
#############################

if !exist('correct_spikes_threshold','var')
    correct_spikes_threshold = [];
end 
if !exist('ccd_filter_para','var')
    ccd_filter_para = [];
end
if !exist('do_ccd_filter','var')
    do_ccd_filter = [];
end

if isempty(ccd_filter_para)
    if isempty(correct_spikes_threshold)
        do_ccd_filter_preference = 'NO';
        ccd_filter_para = sprintf('{"threshold": %f }',correct_spikes_threshold_default);
    elseif isinf(correct_spikes_threshold)
        do_ccd_filter = 'NO';
        ccd_filter_para = sprintf('{"threshold": %f }',correct_spikes_threshold_default);
    else
        do_ccd_filter_preference = 'YES';
        ccd_filter_para = sprintf('{"threshold": %f }',correct_spikes_threshold);
    end
else
    do_ccd_filter_preference = 'YES';
end
if isempty(do_ccd_filter)
    do_ccd_filter = do_ccd_filter_preference;
end

##########################
##### create parfile #####
##########################

fid=fopen(parfilename,'wb');

if (fid==-1)
	disp('problems writing hst-parameterfile ! Check write permissions');
	return
end

# write out parameter file

fprintf(fid,'# HST_SLAVE PARAMETER FILE\n\n');
fprintf(fid,'# Parameters defining the projection file series\n');
fprintf(fid,'FILE_PREFIX = %s/%s\n',direc,prefix);
fprintf(fid,'NUM_FIRST_IMAGE = %d # No. of first projection file\n',num_first_image);
fprintf(fid,'NUM_LAST_IMAGE = %d # No. of last projection file\n',num_last_image);
fprintf(fid,'NUMBER_LENGTH_VARIES = %s\n',number_length_varies);
fprintf(fid,'LENGTH_OF_NUMERICAL_PART = %d # No. of characters\n',length_of_numerical_part);
fprintf(fid,'FILE_POSTFIX = %s\n',file_postfix);
fprintf(fid,'FILE_INTERVAL = %d # Interval between input files\n',file_interval);

fprintf(fid,'\n# Parameters defining the projection file format\n');
fprintf(fid,'NUM_IMAGE_1 = %d # Number of pixels horizontally\n',num_image_1);
fprintf(fid,'NUM_IMAGE_2 = %d # Number of pixels vertically\n',num_image_2);
fprintf(fid,'IMAGE_PIXEL_SIZE_1 = %f # Pixel size horizontally (microns)\n',image_pixel_size_1);
fprintf(fid,'IMAGE_PIXEL_SIZE_2 = %f # Pixel size vertically\n',image_pixel_size_2);

fprintf(fid,'\n# Parameters defining background treatment\n');

if correct
	fprintf(fid,'SUBTRACT_BACKGROUND = %s # Subtract background from data\n',subtract_background);
	fprintf(fid,'BACKGROUND_FILE = %s\n',background_file);
else
	fprintf(fid,'SUBTRACT_BACKGROUND = %s # No background subtraction\n',subtract_background);
	fprintf(fid,'BACKGROUND_FILE = %s\n',not_available);
end
	                               
fprintf(fid,'\n# Parameters defining flat-field treatment\n');
if correct
	fprintf(fid,'CORRECT_FLATFIELD = %s # Divide by flat-field image\n',correct_flatfield);
	fprintf(fid,'FLATFIELD_CHANGING = %s # Series of flat-field files\n',flatfield_changing);
	fprintf(fid,'FLATFIELD_FILE = %s\n',flatfield_file);
	fprintf(fid,'FF_PREFIX = %s\n',ff_prefix);
	fprintf(fid,'FF_NUM_FIRST_IMAGE = %d # No. of first flat-field file\n',ff_num_first_image);
	fprintf(fid,'FF_NUM_LAST_IMAGE = %d # No. of last flat-field file\n',ff_num_last_image);
	fprintf(fid,'FF_NUMBER_LENGTH_VARIES = %s\n',number_length_varies);
	fprintf(fid,'FF_LENGTH_OF_NUMERICAL_PART = %d # No. of characters\n',length_of_numerical_part);
	fprintf(fid,'FF_POSTFIX = %s\n',file_postfix);
	fprintf(fid,'FF_FILE_INTERVAL = %d # Interval between flat-field files\n',ff_file_interval);
else
	fprintf(fid,'CORRECT_FLATFIELD = %s # No flat-field correction\n',correct_flatfield);
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

fprintf(fid,'\nTAKE_LOGARITHM = %s # Take log of projection values\n',take_logarithm);
if !isempty(doubleffcorrection)
    fprintf(fid,'DOUBLEFFCORRECTION = "filtered_mean.edf" # double flatfield\n');
endif

fprintf(fid,'\n# Parameters defining experiment\n');
fprintf(fid,'ANGLE_BETWEEN_PROJECTIONS = %f # Increment angle in degrees\n',angle_between_projections);
fprintf(fid,'ROTATION_VERTICAL = %s\n',rotation_vertical);
fprintf(fid,'ROTATION_AXIS_POSITION = %f # Position in pixels\n',rotation_axis_position);

fprintf(fid,'\n# Parameters defining reconstruction\n');
fprintf(fid,'OUTPUT_SINOGRAMS = %s # Output sinograms to files or not\n',output_sinograms);
fprintf(fid,'OUTPUT_RECONSTRUCTION = %s # Reconstruct and save or not\n',output_reconstruction);
fprintf(fid,'START_VOXEL_1 =   %4d # X-start of reconstruction volume\n',start_voxel_1);
fprintf(fid,'START_VOXEL_2 =   %4d # Y-start of reconstruction volume\n',start_voxel_2);
fprintf(fid,'START_VOXEL_3 =   %4d # Z-start of reconstruction volume\n',start_voxel_3);
fprintf(fid,'END_VOXEL_1 =   %4d # X-end of reconstruction volume\n',end_voxel_1);
fprintf(fid,'END_VOXEL_2 =   %4d # Y-end of reconstruction volume\n',end_voxel_2);
fprintf(fid,'END_VOXEL_3 =   %4d # Z-end of reconstruction volume\n',end_voxel_3);
fprintf(fid,'OVERSAMPLING_FACTOR = %d # 0 = Linear, 1 = Nearest pixel\n',oversampling_factor);
fprintf(fid,'ANGLE_OFFSET = %f # Reconstruction rotation offset angle in degrees\n',angle_offset);
fprintf(fid,'CACHE_KILOBYTES = %d # Size of processor cache (L2) per processor (KBytes)\n',cache_kilobytes);
fprintf(fid,'SINOGRAM_MEGABYTES = %d # Maximum size of sinogram storage (megabytes)\n',sinogram_megabytes);

# adding extra features for PyHST
fprintf(fid,'\n# Parameters extra features PyHST\n');
fprintf(fid,'DO_CCD_FILTER = %s # CCD filter (spikes)\n',do_ccd_filter);
fprintf(fid,'CCD_FILTER = "%s"\n',ccd_filter);
fprintf(fid,'CCD_FILTER_PARA = %s\n',ccd_filter_para);
fprintf(fid,'DO_SINO_FILTER = %s # Sinogram filter (rings)\n',do_sino_filter);
fprintf(fid,'SINO_FILTER = "%s"\n','SINO_Filter');
fprintf(fid,'ar = Numeric.ones(%d,''f'')\n',num_image_1);
fprintf(fid,'ar[0]=0.0\n');
if (correct_rings_nb > 0)
	fprintf(fid,'ar[2:%d]=0.0\n',2*(correct_rings_nb+1));
end
fprintf(fid,'SINO_FILTER_PARA = {"FILTER": ar }\n');
fprintf(fid,'DO_AXIS_CORRECTION = %s # Axis correction\n',do_axis_correction);
fprintf(fid,'AXIS_CORRECTION_FILE = %s\n',axis_correction_file);
fprintf(fid,'OPTIONS= { ''padding'':''%s'' , ''axis_to_the_center'':''%s''} # Padding and position axis\n',padding,axis_to_the_center);   
if !isempty(zeroclipvalue)
    fprintf(fid,'ZEROCLIPVALUE = %g # Minimum value of radiographs after flat / before log\n',zeroclipvalue);    
endif
if !isempty(oneclipvalue)
    fprintf(fid,'ONECLIPVALUE = %g # Maximum value of radiographs after flat / before log\n',oneclipvalue);    
endif
if !strcmp(do_projection_median,'NO')
    fprintf(fid,'DO_PROJECTION_MEDIAN = %s # Calculate median of all projections\n',do_projection_median);
    fprintf(fid,'PROJECTION_MEDIAN_FILENAME = "%s" # Name output file median calculation\n',projection_median_filename)
endif
if !strcmp(do_projection_mean,'NO')
    fprintf(fid,'DO_PROJECTION_MEAN = %s # Calculate median of all projections\n',do_projection_mean);
    fprintf(fid,'PROJECTION_MEDIAN_FILENAME = "%s" # Name output file mean calculation\n',projection_mean_filename)
endif

fprintf(fid,'\n# Parameters defining output file / format\n');
fprintf(fid,'OUTPUT_FILE = %s\n',output_file);

fprintf(fid,'\n# Reconstruction program options\n');
fprintf(fid,'DISPLAY_GRAPHICS = %s # No images\n',display_graphics);

fclose(fid);

if (dir(parfilename).statinfo.uid == getuid)
    # we change permissions if the file belongs to us
    system(['chmod 666 ' parfilename]);  
endif












