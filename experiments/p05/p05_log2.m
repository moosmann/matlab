function [p,cur] = p05_log( file )
% Read log file of beamline P05.
%
% ARGUMENTS
% file : string. Folder of or path to log file
% 
% Written by Julian Moosmann, 2016-12-12. Last version: 2016-12-12
%
% p = p05_log(filename)

%% Default arguments %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 1
    file = '/asap3/petra3/gpfs/p05/2015/data/11001102/raw/hzg_wzb_mgag_14/hzg_wzb_mgag_14scan.log';
    %file = '/asap3/petra3/gpfs/p05/2016/data/11001978/raw/mah_28_15R_top/mah_28_15R_topscan.log';
end

%% Main %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if exist( file , 'file') ~= 2
    s = dir( [file '/*p.log'] );    
    file = [s.folder filesep s.name];
end
fid = fopen( file );

%% Read log file until ring current
cl = {{[]}};
cl_counter = 1;
while 1
    c = textscan( fid, '%s', 1, 'Delimiter', {'\n', '\r'} );
    s = c{1}{1};    
    if ~strcmpi( s(1), '*')
        cl{1}{cl_counter} = s;
        cl_counter = cl_counter + 1;
    end    
    if strcmpi( c{1}{1}(1:6), '/PETRA' )
        c_break = c;
        break
    end    
end

%% Possible parameter to read out

% p info
p.Energy = '%f';
p.undulator_gap = '%f';
p.undulator_harmonic = '%u';
p.n_dark = '%u';
p.n_ref_min = '%u';
p.n_ref_max = '%u';
p.n_img = '%u';
p.n_angle = '%u';
p.ref_count = '%u';
p.img_bin = '%u';
p.img_xmin = '%u';
p.img_xmax = '%u';
p.img_ymin = '%u';
p.img_ymax = '%u';
p.exptime = '%u';

% camera info 
p.ccd_pixsize = '%f';
p.ccd_xsize =  '%u';
p.ccd_ysize =  '%u';
p.o_screen_changer = '%f';
p.o_focus = '%f';
p.o_aperture = '%f';
p.o_lens_changer = '%f';
p.o_ccd_high = '%f';
p.magn = '%f';

% apparatus info
p.pos_s_pos_x = '%f';
p.pos_s_pos_y = '%f';
p.pos_s_pos_z = '%f';
p.pos_s_stage_z = '%f';
p.s_in_pos = '%f';
p.s_out_dist = '%f';
p.o_ccd_dist = '%f';

% optimization info
p.p05_dcm_xtal2_pitch_delta = '%f';
p.com_delta_threshhold =  '%f';

%% Compare possible paramters with strings from log file
fn = fieldnames( p );
for nn = 1:numel( fn )
    p_field = fn{nn};        
    b = contains( cl{1}, fn(nn) );
    if sum( b ) == 1
        p_val = p.(p_field);
        c = textscan( cl{1}{b}, sprintf( '%s %s', p_field, p_val ), 1, 'Delimiter', {'=', ':'});  
        p.(p_field) = c{1};
    elseif sum( b ) > 1
        fprintf( 'WARNING: Parameter found more than once' );
    else
        p = rmfield( p, p_field);        
    end
    
end

%% ring current

% ************************************************
% /PETRA/Idc/Buffer-0/I.SCH
% @1480217687541[99.99254@1480217687487]
% dark /gpfs/current/raw/mah_28_15R_top/mah_28_15R_top00000.dar
% /PETRA/Idc/Buffer-0/I.SCH
% @1480217687953[99.98162@1480217687928]    

% ************************************************
% /PETRA/Idc/Buffer-0/I.SCH
% @1436994685507[89.40661@1436994685435]
% 
% ref /gpfs/current/raw/hzg_wzb_mgag_00/hzg_wzb_mgag_0000005.ref       0.00000
% /PETRA/Idc/Buffer-0/I.SCH
% @1436994686180[89.39424@1436994686098]

% Move back to first occurence of /PETRA...
fseek( fid, - ( length( c_break{1}{1} ) + 1 ), 'cof' );

% Read all current fields
fs = '/%*s@%*f[%f@%*f]%s/%*s@%*f[%f@%*f]%*s';
c  = textscan( fid, fs, 'Delimiter', {'\n', '\r'},  'MultipleDelimsAsOne', 1);

num_files = numel( c{1} );
num_proj = p.n_angle + 1;
num_ref = num_files - ( num_proj + p.n_dark );
ref_count = 1;
proj_count = 1;
for nn = 1:num_files
    str = c{2}{nn};
    val = ( c{1}(nn) + c{3}(nn) ) / 2;    
     switch str(1:3)
         case 'ref' 
             st_c = textscan(str, '%*s%s%f');
             cur.ref(ref_count).ind = str2double(st_c{1}{1}(end-8:end-4));
             cur.ref(ref_count).angle = st_c{2};
             cur.ref(ref_count).val = val;
             ref_count = ref_count + 1;
         case 'img'
             st_c = textscan(str, '%*s%s%f');
             cur.proj(proj_count).ind = str2double(st_c{1}{1}(end-8:end-4));
             cur.proj(proj_count).angle = st_c{2};
             cur.proj(proj_count).val = val;
             proj_count = proj_count + 1;
     end
end

fclose( fid );
% End %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%