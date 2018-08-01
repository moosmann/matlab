function [par, cur, cam] = p05_log( file )
% Read log file of beamline P05.
%
% ARGUMENTS
% file : string. Folder of or path to log file
%
% Written by Julian Moosmann, 2016-12-12. Last version: 2016-12-12
%
% [par, cur, cam] = p05_log(filename)

%% Default arguments %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 1
    file = '/asap3/petra3/gpfs/p05/2016/data/11001994/raw/szeb_74_13/scan.log';
end

%% Main %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if exist( file , 'file') ~= 2
    s = dir( [file '/*scan.log'] );
    file = [s.folder filesep s.name];
end

[~, name] = fileparts( file );
if length( name ) > 4
    cam = 'EHD';
    [par, cur] = EHD_log( file );
    if isfield( par, 'ccd_xsize') && isequal( par.ccd_xsize, 5120 )
        cam = 'KIT';
    end
else
    cam = 'KIT';
    [par, cur] = KIT_log( file );
end

%% Fix inconsistent naming schemes
if isfield( par, 'Energy')
    par.energy = par.Energy;
    par = rmfield( par, 'Energy');
end

if isfield( par, 'eff_pix_size')
    par.eff_pixel_size = abs( par.eff_pix_size * 1e-3 );
    par = rmfield( par, 'eff_pix_size');
elseif isfield( par, 'eff_pix')
    par.eff_pixel_size = abs( par.eff_pix * 1e-3 );
    par = rmfield( par, 'eff_pix');
elseif isfield( par, 'ccd_pixsize') && isfield( par, 'magn')
    par.eff_pixel_size = abs( par.ccd_pixsize / par.magn * 1e-3 );
    par = rmfield( par, 'ccd_pixsize');
end

if isfield( par, 'camera_distance')
    par.sample_detector_distance = par.camera_distance / 1000;
        par = rmfield( par, 'camera_distance');
elseif isfield( par, 'camera_dist')
    par.sample_detector_distance = par.camera_dist / 1000;
    par = rmfield( par, 'camera_dist');
elseif isfield( par, 'o_ccd_dist')
    par.sample_detector_distance = par.o_ccd_dist / 1000;
    par = rmfield( par, 'o_ccd_dist');
end

if isfield( par, 'num_projections' )
    par.num_proj = double( par.num_projections );
elseif isfield( par, 'projections' )
    par.num_proj = double( par.projections );
elseif isfield( par, 'n_angles' )
    par.num_proj = double( par.n_angles );
elseif isfield( par, 'n_angle' )
    par.num_proj = double( par.n_angle );
end

if isfield( par, 'exptime' )
    par.exposure_time = double( par.exptime );
    par = rmfield( par, 'exptime' );
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [par, cur] = KIT_log( file )
fid = fopen( file );
c = textscan( fid, '%s');
fclose( fid );

%% Possible parameter list
par.energy = '%f';
par.num_flat_per_block = '%u' ;
par.ROI = '%u%u%u%u';
par.scan_duration = '%f';
par.eff_pix_size = '%f';
par.eff_pix = '%f';
par.projections = '%u';
par.num_projections = '%u';
par.ref_prefix = '%s';
par.sample = '%s';
par.rotation = '%f';
par.num_img_per_block = '%u';
par.angle_list = '%f';
par.angle_order = '%s';
par.height_steps = '%u';
par.dark_prefix = '%s';
par.num_dark_img = '%u';
par.pos_dpc_s_pos_z = '%f';
par.camera_distance = '%f';
par.camera_dist = '%f';
par.proj_prefix = '%s';
par.exposure_time = '%f';
par.off_axes = '%u';

%% Compare possible paramters with strings from log file
fn = fieldnames( par );
for nn = 1:numel( fn )
    p_field = fn{nn}; % cell array
    b = contains( c{1}, fn(nn) ); %
    if sum( b ) == 1
        p_val = par.(p_field);
        %         cl = textscan( c{1}{b}, sprintf( '%s %s', p_field, p_val ), ...
        %             1, 'Delimiter', {'=',','}, 'CollectOutput', 1, 'MultipleDelimsAsOne', 1);
        cl = textscan( c{1}{b}(2+numel(p_field):end), p_val, ...
            'Delimiter', {'=',','}, 'CollectOutput', 1, 'MultipleDelimsAsOne', 1);
        
        val = cl{1};
        if iscell( val )
            par.(p_field) = val{1};
        else
            par.(p_field) = val;
        end
    elseif sum( b ) > 1
        fprintf( 'WARNING: Parameter not unique' );
    else
        par = rmfield( par, p_field);
    end
    
end

%% Ring current
[folder] = fileparts( file );
proj_cur = single( read_dat( sprintf('%s/proj_current.dat', folder ) ) );
ref_cur = single( read_dat( sprintf('%s/ref_current.dat', folder ) ) );

if ~isempty(par.projections)
    num_proj = single(par.projections);
elseif ~isempty(par.num_projections)
    num_proj = single(par.num_projections);
end
if ~isfield( par, 'angle_list' )
    par.angle_list = 180 * (0:num_proj - 1) / num_proj;
end
for nn = 1:numel( proj_cur )
    cur.proj(nn).val = proj_cur(nn);
    cur.proj(nn).ind = nn - 1;
    cur.proj(nn).angle = par.angle_list(nn);
end
for nn = 1:numel( ref_cur )
    cur.ref(nn).val = ref_cur(nn);
    cur.ref(nn).ind = nn - 1;
    cur.ref(nn).angle = 0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [par, cur] = EHD_log( file )
fid = fopen( file );
%% Read log file until ring current
cl = {{[]}};
cl_counter = 1;
while 1
    c = textscan( fid, '%s', 1, 'Delimiter', {'\n', '\r'} );
    if ~isempty( c{1} )
        s = c{1}{1};
        if ~isempty( s ) && ~strcmpi( s(1), '*')
            cl{1}{cl_counter} = s;
            cl_counter = cl_counter + 1;
        end
        if ~isempty( c{1}{1} ) && strcmpi( c{1}{1}(1:6), '/PETRA' )
            c_break = c;
            break
        end
    else
        c_break = [];
        break
    end
end

%% Possible parameters to read

% setup info
par.Energy = '%f';
par.undulator_gap = '%f';
par.undulator_harmonic = '%u';
par.n_dark = '%u';
par.n_ref_min = '%u';
par.n_ref_max = '%u';
par.n_img = '%u';
par.n_angle = '%u';
par.ref_count = '%u';
par.img_bin = '%u';
par.img_xmin = '%u';
par.img_xmax = '%u';
par.img_ymin = '%u';
par.img_ymax = '%u';
par.exptime = '%u';

% camera info
par.ccd_pixsize = '%f';
par.ccd_xsize =  '%u';
par.ccd_ysize =  '%u';
par.o_screen_changer = '%f';
par.o_focus = '%f';
par.o_aperture = '%f';
par.o_lens_changer = '%f';
par.o_ccd_high = '%f';
par.magn = '%f';
par.eff_pix = '%f';

% apparatus info
par.pos_s_pos_x = '%f';
par.pos_s_pos_y = '%f';
par.pos_s_pos_z = '%f';
par.pos_s_stage_z = '%f';
par.s_in_pos = '%f';
par.s_out_dist = '%f';
par.o_ccd_dist = '%f';

% optimization info
par.p05_dcm_xtal2_pitch_delta = '%f';
par.com_delta_threshhold =  '%f';

%% Compare possible paramters with strings from log file
fn = fieldnames( par );
for nn = 1:numel( fn )
    p_field = fn{nn};
    b = contains( cl{1}, fn(nn) );
    if sum( b ) == 1
        p_val = par.(p_field);
        c = textscan( cl{1}{b}, sprintf( '%s %s', p_field, p_val ), 1, 'Delimiter', {'=', ':'});
        par.(p_field) = c{1};
    elseif sum( b ) > 1
        fprintf( 'WARNING: Parameter found more than once' );
    else
        par = rmfield( par, p_field);
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
if ~isempty( c_break )
    fseek( fid, - ( length( c_break{1}{1} ) + 1 ), 'cof' );
    
    % Read all current fields
    fs = '/%*s@%*f[%f@%*f]%s/%*s@%*f[%f@%*f]%*s';
    c  = textscan( fid, fs, 'Delimiter', {'\n', '\r'},  'MultipleDelimsAsOne', 1);
    
    ref_count = 1;
    proj_count = 1;
    for nn = 1:numel( c{1} )
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
else
    cur = [];
end

fclose( fid );
% End %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%