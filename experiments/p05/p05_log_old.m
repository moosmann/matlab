function [scan, cam, app, opt, ref, proj] = p05_log_old( file )
% Read log file of beamline P05.
%
% ARGUMENTS
% file : string. Folder of or path to log file
% 
% Written by Julian Moosmann, 2016-12-12. Last version: 2016-12-12
%
% [scan, cam, app, opt, cur] = p05_log(filename)

%% Default arguments %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 1
    %file = '/asap3/petra3/gpfs/p05/2015/data/11001102/raw/hzg_wzb_mgag_14/hzg_wzb_mgag_14scan.log';
    %file = '/asap3/petra3/gpfs/p05/2015/data/11001102/raw/hzg_wzb_mgag_14/';
    file = '/asap3/petra3/gpfs/p05/2016/data/11001978/raw/mah_28_15R_top/mah_28_15R_topscan.log';
end

%% Main %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if exist( file , 'file') ~= 2
    s = dir( [file '/*scan.log'] );    
    file = [s.folder filesep s.name];
end
fid = fopen( file );

%% scan info
textscan( fid, '%s', 2, 'Delimiter', '\n' ); % Skip first two lines
scan.Energy = '%f';
scan.undulator_gap = '%f';
scan.undulator_harmonic = '%u';
scan.n_dark = '%u';
scan.n_ref_min = '%u';
scan.n_ref_max = '%u';
scan.n_img = '%u';
scan.n_angle = '%u';
scan.ref_count = '%u';
scan.img_bin = '%u';
scan.img_xmin = '%u';
scan.img_xmax = '%u';
scan.img_ymin = '%u';
scan.img_ymax = '%u';
scan.exptime = '%u';

fn = fieldnames( scan );
for nn = 1:numel( fn )
    s = fn{nn};
    v = scan.(s);   
    c = textscan( fid, sprintf( '%s=%s', s, v ) );
    scan.(s) = c{1};
    if isempty( scan.(s) )
        scan = rmfield( scan, s);
    end
end
textscan( fid, '%s', 1, 'Delimiter', '\n' );
textscan( fid, '%s', 1, 'Delimiter', '\n' );

%% camera info 
cam.ccd_pixsize = '%f';
cam.ccd_xsize =  '%u';
cam.ccd_ysize =  '%u';
cam.o_screen_changer = '%f';
cam.o_focus = '%f';
cam.o_aperture = '%f';
cam.o_lens_changer = '%f';
cam.o_ccd_high = '%f';
cam.magn = '%f';

fn = fieldnames( cam );
for nn = 1:numel( fn )
    s = fn{nn};
    v = cam.(s);    
    c = textscan( fid, sprintf( '%s:%s', s, v ) );
    cam.(s) = c{1};
end

textscan( fid, '%s', 4, 'Delimiter', '\n');
textscan( fid, '%s', 2, 'Delimiter', '\n');

%% apparatus info
app.pos_s_pos_x = '%f';
app.pos_s_pos_y = '%f';
app.pos_s_pos_z = '%f';
app.pos_s_stage_z = '%f';
app.s_in_pos = '%f';
app.s_out_dist = '%f';
app.o_ccd_dist = '%f';

fn = fieldnames( app );
for nn = 1:numel( fn )
    s = fn{nn};
    v = app.(s);
    % fprintf( 'fieldname: %s : ', s ); disp( v )   
    c = textscan( fid, sprintf( '%s=%s', s, v ), 'Delimiter', '\n' );
    app.(s) = c{1};    
    if isempty( app.(s) )
        app = rmfield( app, s);
    end
end
textscan( fid, '%s', 1, 'Delimiter', '\n');

%% optimization info
opt.p05_dcm_xtal2_pitch_delta = '%f';
opt.com_delta_threshhold =  '%f';

fn = fieldnames( opt );
for nn = 1:numel( fn )
    s = fn{nn};
    v = opt.(s);
    c = textscan( fid, sprintf( '%s=%s', s, v ), 'Delimiter', '\n' );
    opt.(s) = c{1};
end

%% ring current
if isfield(scan, 'undulator_gap')
% ************************************************
% /PETRA/Idc/Buffer-0/I.SCH
% @1480217687541[99.99254@1480217687487]
% dark /gpfs/current/raw/mah_28_15R_top/mah_28_15R_top00000.dar
% /PETRA/Idc/Buffer-0/I.SCH
% @1480217687953[99.98162@1480217687928]    
    formatSpec = '%*s\n%*s\n@%*f[%f@%f]\n%s /gpfs/current/raw/%s %f\n%*s\n@%*f[%f@%f]%*s\n';
else
% ************************************************
% /PETRA/Idc/Buffer-0/I.SCH
% @1436994685507[89.40661@1436994685435]
% 
% ref /gpfs/current/raw/hzg_wzb_mgag_00/hzg_wzb_mgag_0000005.ref       0.00000
% /PETRA/Idc/Buffer-0/I.SCH
% @1436994686180[89.39424@1436994686098]
    formatSpec = '%*s\n%*s\n@%*f[%f@%f]\n\n%s /gpfs/current/raw/%s %f\n%*s\n@%*f[%f@%f]\n%*s\n';    
end
c  = textscan( fid, formatSpec);
num_files = numel( c{1} );
num_proj = scan.n_angle + 1;
num_ref = num_files - ( num_proj + scan.n_dark );
ref(num_ref).name = '';
ref_count = 1;
proj(num_proj).name = '';
proj_count = 1;
for nn = 1:num_files
    [~, name, ~] = fileparts( c{4}{nn} );
    ind = str2double( name(end-4:end) );
    nnangle = c{5}(nn);
    cur = ( c{1}(nn) + c{6}(nn) ) / 2;
     switch c{3}{nn} 
         case 'ref'            
             ref(ref_count).name = name;
             ref(ref_count).ind = ind;
             ref(ref_count).angle = nnangle;
             ref(ref_count).cur = cur;
             ref_count = ref_count + 1;
         case 'img'
             proj(proj_count).name = name;
             proj(proj_count).ind = ind;
             proj(proj_count).angle = nnangle;
             proj(proj_count).cur = cur;
             proj_count = proj_count + 1;
     end
end

% End of Scan
% End of File

fclose( fid );
% End %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%