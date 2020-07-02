
%% Parameters
clear all
read_scans = 1;
ignore_case = 1;
process_tags = {'cpd'};
stitch_only = 0;
run_height_scans = 1;
run_single_height_scans = 1;
% Empiric values from full range histograms and median filtering
bin_limits = [-0.0018 0.0092]; 
% Output directory
outpath = '/asap3/petra3/gpfs/p05/2016/data/11001978/processed/segmentation/';
%outpath = '/asap3/petra3/gpfs/p05/2017/data/11003288/processed/segmentation/';

%% Beamtimes
beamtime = {
    2016 11001978
    2017 11003950
    2017 11003288
    2017 11003440
    2017 11003773
    2017 11004016
    2018 11004263
    2018 11004936
    2018 11005553
    2019 11005842
    2019 11006704
    2019 11006991
    };

samples = synchroload_samples;
all_counter = 0;

%% Loop over cpd and load scans
for exp_type = { 'cpd', 'load'}
    et = exp_type{1};
    et_counter = 0;
    
    %% Loop over implant types
    for imp_type = { 'Mg5Gd', 'Mg10Gd', 'Ti', 'PEEK' }
        it = imp_type{1};
        it_counter = 0;
        fprintf( '\n\n%s %s', et, it )
        name_pattern = samples.(et).(it);
        % add underscore for unique match using regexp
        for np = 1:numel( name_pattern )
            name_pattern{np} = [ '_' name_pattern{np}];
            sample_name = [ 'c' name_pattern{np} ];
            samples.control.(et).(it).( sample_name ) = 0;
        end
        
        %% Loop over beamtimes
        for bt = 1:size( beamtime, 1 )
            % beamtime folder
            if size( beamtime, 1 ) == 1
                beamtime_year = beamtime{1};
                beamtime_id = beamtime{2};
            else
                beamtime_year = beamtime{bt,1};
                beamtime_id = beamtime{bt,2};
            end
            proc_path = sprintf( '/asap3/petra3/gpfs/p05/%u/data/%u/processed', beamtime_year, beamtime_id );
            % scan folder
            dir_struct = dir( proc_path );
            match2 = contains( {dir_struct.name}, name_pattern, 'IgnoreCase', ignore_case );
            match = zeros( [1, numel(dir_struct)], 'logical' );
            sample_arr = cell( [1, numel(dir_struct)] );
            for mm = 1:numel( match )
                for np = 1:numel(name_pattern )
                    name_pattern_match = contains( dir_struct(mm).name, name_pattern{np}, 'IgnoreCase', ignore_case );
                    sample_name = [ 'c' name_pattern{np} ];
                    if name_pattern_match
                        match(mm) =  1;
                        samples.control.(et).(it).( sample_name ) = 1;
                        sample_arr{mm} = name_pattern{np};
                    end
                end
            end
            if ~isempty( match ) && ~isempty( match2 ) && ~isequal( match, match2 )
                error( '\n Something is wrong' )
            end
            scan_struct = dir_struct(match);
            fprintf( '\nexp %2u %s %s %u %u %u: %s found', bt, et, it, beamtime_year, beamtime_id, length( scan_struct ) );
            sample_arr = sample_arr(match);
            
            % Loop over scans
            for mm = 1:length( scan_struct )
                folder_name = scan_struct(mm).name;
                scan_path = [ scan_struct(mm).folder filesep folder_name ];
                sample = sample_arr{mm};
                % Reset parameters
                data_path = [];
                reco_log = [];
                effective_pixel_size = [];
                effective_pixel_size_binned = [];
                energy = [];
                camera = [];
                im_shape_raw = [];
                bin = [];
                raw_roi = [];
                raw_path = [];
                scan_log = [];
                s_stage_z = [];
                
                if isfolder( scan_path )
                    reco_path = [scan_path filesep 'reco' ];
                    if exist( reco_path, 'dir')
                        %% Find float recos
                        d = dir( [ reco_path '/float*' ] );
                        disdir =  [ d(:).isdir ];
                        floats_namecell = { d( disdir ).name };
                        if isempty( floats_namecell )
                            continue
                        end
                        %% Choose first float reco
                        subfolder_reco = floats_namecell{1};
                        data_path = [reco_path filesep subfolder_reco];
                        
                        %% Read reco log file
                        dir_struct = dir( [reco_path '/reco*.log'] );
                        reco_log = [ reco_path  filesep  dir_struct(1).name ];
                        fid = fopen( reco_log );
                        c = textscan( fid, '%s', 'Delimiter', {'\n', '\r'} );
                        c = c{1};
                        fclose( fid );
                        %% Effective pixelsize
                        for ll = 1:numel( c )
                            t = regexp( c{ll}, 'effective_pixel_size' );
                            if t
                                cc = textscan( c{ll}, '%*s : %f %*s');
                                effective_pixel_size = cc{1};
                                break
                            end
                        end
                        %% Energy
                        for ll = 1:numel( c )
                            t = regexp( c{ll}, 'energy' );
                            if t
                                cc = textscan( c{ll}, '%*s : %f %*s');
                                energy = cc{1};
                                break
                            end
                        end
                        %% Camera
                        for ll = 1:numel( c )
                            t = regexp( c{ll}, 'camera' );
                            if t
                                cc = textscan( c{ll}, '%*s : %s');
                                camera = cc{1}{1};
                                break
                            end
                        end
                        %% im shape raw
                        for ll = 1:numel( c )
                            t = regexp( c{ll}, {'im_shape_raw', 'raw_image_shape', 'raw_image_shape_raw'} );
                            %disp( c{ll} )
                            if sum( [t{:}] )
                                cc = textscan( c{ll}, '%*s : %u %u');
                                im_shape_raw = [cc{:}];
                                break
                            end
                        end
                        %% Binning factor
                        for ll = 1:numel( c )
                            t = regexp( c{ll}, {'raw_bin', 'raw_binning_factor', 'raw_binning'} );
                            if sum( [t{:}] )
                                cc = textscan( c{ll}, '%*s : %u %u');
                                bin = [cc{:}];
                                break
                            end
                        end
                        if ~isscalar( bin )
                            t = regexp( subfolder_reco, 'bin|Bin' );
                            if ~isempty( t )
                                bin = str2double( subfolder_reco(t+3) );
                            else
                                bin = 1;
                            end
                        end
                        %% Raw roi
                        for ll = 1:numel( c )
                            t = regexp( c{ll}, 'raw_roi' );
                            if t
                                cc = textscan( c{ll}, '%*s : %u %u');
                                raw_roi = [ cc{1} cc{2} ];
                                break
                            end
                        end
                        if ~exist( 'raw_roi', 'var' ) || ( exist( 'raw_roi', 'var' ) && isempty( raw_roi ) )
                            raw_roi(1) = 1;
                            raw_roi(2) = im_shape_raw(2);
                        end
                        if numel( raw_roi ) ~=2
                            edit( reco_log )
                            error( '\n raw_roi paramter is wrong' )
                        end
                        %% Raw path
                        for ll = 1:numel( c )
                            t = regexp( c{ll}, 'scan_path' );
                            if t
                                cc = textscan( c{ll}, '%*s : %s');
                                raw_path = cc{1}{1};
                                break
                            end
                        end
                        %% Scan log
                        scan_log = [raw_path filesep folder_name 'scan.log'];
                        if ~exist( scan_log, 'file' )
                            scan_log = [raw_path filesep 'scan.log'];
                        end
                        if ~exist( scan_log, 'file' )
                            scan_log = dir( [raw_path filesep '*scan.log'] );
                            scan_log = [scan_log.folder filesep scan_log.name];
                        end
                        if ~exist( scan_log, 'file' )
                            error( '\nScan log not found!\n' )
                        end
                        fid = fopen( scan_log );
                        cell_of_lines = textscan( fid, '%s', 'Delimiter', {'\n', '\r'} );
                        cell_of_lines = cell_of_lines{1};
                        fclose( fid );
                        %% s_stage_z
                        for ll = 1:numel( cell_of_lines )
                            t = regexp( cell_of_lines{ll}, 's_stage_z' );
                            if t
                                cc = textscan( cell_of_lines{ll}, '%s', 'Delimiter', {'=',',',' '}, 'CollectOutput', 1, 'MultipleDelimsAsOne', 1);
                                s_stage_z = str2double( cc{1}{2} ) * 1000;
                                break
                            end
                        end
                        %% effective_pixel_size binned
                        effective_pixel_size_binned = double( bin ) * effective_pixel_size;
                        %% Unique and descriptive name
                        name = sprintf( '%s_%s_E%05.0feV_cam%ux%u_pix%04.0fnm%s_%u_%u_%s', et, it, energy, im_shape_raw, 1000*effective_pixel_size, sample,  beamtime_year, beamtime_id, folder_name);
                        %% Create parameter struct
                        tmp.type = et;
                        tmp.material = it;
                        tmp.energy = energy;
                        tmp.im_shape_raw = im_shape_raw;
                        tmp.raw_roi = raw_roi;
                        tmp.effective_pixel_size = effective_pixel_size;
                        tmp.effective_pixel_size_binned = effective_pixel_size_binned;
                        tmp.bin = bin;
                        tmp.sample = sample(2:end);
                        tmp.beamtime_year = beamtime_year;
                        tmp.beamtime_id = beamtime_id;
                        tmp.folder_name = folder_name;
                        tmp.camera = camera;
                        tmp.name = name;
                        tmp.s_stage_z = s_stage_z;
                        tmp.raw_path = raw_path;
                        tmp.data_path = data_path;
                        tmp.proc_path = proc_path;
                        tmp.reco_log = reco_log;
                        tmp.scan_log = scan_log;
                        
                        % Assign to groups
                        it_counter = it_counter + 1;
                        scans.(et).(it)(it_counter) = tmp;
                        et_counter = et_counter + 1;
                        scans.(et).all(et_counter) = tmp;
                        all_counter = all_counter + 1;
                        scans.all(all_counter) = tmp;
                        fprintf( '\n %2u. bin: %u, effpix: %f, cam: %s, shape: %u %u, name: %s', mm, bin, effective_pixel_size, camera, im_shape_raw, name )
                    end
                end
            end
        end
    end % imp_type
end % exp_type
%% Save parameter struct
CheckAndMakePath( outpath )
mfile_path = [outpath 'matlab/'];
CheckAndMakePath( mfile_path )
filename = [ mfile_path 'scans.mat'];
save( filename, 'scans' );
filename = [ mfile_path 'samples.mat'];
save( filename, 'samples' );
fprintf( '\nFinished reading scans.\n' )

%         csvfilename=[out_path filesep struct_scans(nn).name '_waiting.csv'];
%         cHeader = {'wait_time' 'wait_force' 'wait_pusher'}; %dummy header
%         commaHeader = [cHeader;repmat({','},1,numel(cHeader))]; %insert commaas
%         commaHeader = commaHeader(:)';
%         textHeader = cell2mat(commaHeader); %cHeader in text with commas
%         %write header to file
%         fid = fopen(csvfilename,'w');
%         fprintf(fid,'%s\n',textHeader);
%         fclose(fid);
%         %write data to end of file
%         dlmwrite(csvfilename,datvalues,'-append');

%% Data sets to process %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf( '\n\nData sets to process' )
% Loop over scan types
for nn = 1:numel( process_tags )
    tag = process_tags{nn};
    fprintf( '\n Scan type: %s', tag );
    scans_proc = scans.(tag).all;
    % Loop over scans
    for mm = 1:numel( scans_proc )
        s = scans_proc(mm);
        scan_name = s.name;
        fprintf( '\n %s', s.name )
    end
end
overviewpath = [outpath 'overview'];
CheckAndMakePath( overviewpath );
parpath = [outpath 'scans/'];
CheckAndMakePath( parpath );

%% List : ALL CPD SCANS
fprintf( '\n\n LISTS' )
cpd_list = {
    'cpd_Mg5Gd_E34014eV_cam3056x1621_pix2427nm_77L_2017_11003950_syn22_77L_Mg5Gd_8w_a' 1
    'cpd_Mg5Gd_E34014eV_cam3056x1621_pix2427nm_77L_2017_11003950_syn22_77L_Mg5Gd_8w_b' 1
    'cpd_Mg5Gd_E34014eV_cam3056x1621_pix2427nm_80L_2017_11003950_syn22_80L_Mg5Gd_8w_a' 2
    'cpd_Mg5Gd_E34014eV_cam3056x1621_pix2427nm_80L_2017_11003950_syn22_80L_Mg5Gd_8w_b' 2
    'cpd_Mg5Gd_E34014eV_cam3056x1621_pix2427nm_88R_2017_11003950_syn22_88R_Mg5Gd_4w_a' 3
    'cpd_Mg5Gd_E34014eV_cam3056x1621_pix2427nm_88R_2017_11003950_syn22_88R_Mg5Gd_4w_b' 3
    'cpd_Mg5Gd_E34014eV_cam3056x1621_pix2427nm_99L_2017_11003950_syn22_99L_Mg5Gd_4w_a' 4
    'cpd_Mg5Gd_E34014eV_cam3056x1621_pix2427nm_99L_2017_11003950_syn22_99L_Mg5Gd_4w_b' 4
    'cpd_Mg5Gd_E38374eV_cam3056x3056_pix2403nm_62L_2017_11003440_syn22_62L_Mg5Gd_12w' 0
    'cpd_Mg5Gd_E38374eV_cam3056x3056_pix2403nm_63L_2017_11003440_syn23_63L_Mg5Gd_12w' 0
    'cpd_Mg5Gd_E38374eV_cam3056x3056_pix2403nm_78L_2017_11003440_syn24_78L_Mg5Gd_8w' 0
    'cpd_Mg5Gd_E38402eV_cam3056x3056_pix2403nm_78L_2017_11003440_syn29_78L_Mg5Gd_8w' 0
    'cpd_Mg5Gd_E38402eV_cam3056x3056_pix2403nm_63L_2017_11003440_syn36_63L_Mg5Gd_12w' 0
    'cpd_Mg5Gd_E38402eV_cam3056x3056_pix2403nm_69L_2017_11003440_syn37_69L_Mg10Gd_12w' 0
    'cpd_Mg5Gd_E38402eV_cam3056x3056_pix2403nm_75L_2017_11003440_syn39_75L_Mg5Gd_8w' 0
    'cpd_Mg5Gd_E38402eV_cam3056x3056_pix2403nm_69L_2017_11003440_syn40_69L_Mg10Gd_12w' 0
    'cpd_Mg5Gd_E38402eV_cam3056x3056_pix2403nm_63L_2017_11003440_syn41_63L_Mg5Gd_12w' 0
    'cpd_Mg5Gd_E38402eV_cam3056x3056_pix2403nm_66L_2017_11003440_syn44_66L_Mg5Gd_12w' 0
    'cpd_Mg5Gd_E38402eV_cam3056x3056_pix2403nm_88R_2017_11003440_syn46_88R_Mg5Gd_4w' 0
    'cpd_Mg5Gd_E38402eV_cam3056x3056_pix2403nm_89L_2017_11003440_syn48_89L_Mg5Gd_4w' 0
    'cpd_Mg5Gd_E38402eV_cam3056x3056_pix2403nm_80L_2017_11003440_syn49_80L_Mg5Gd_8w' 0
    'cpd_Mg5Gd_E38402eV_cam3056x3056_pix2403nm_99L_2017_11003440_syn50_99L_Mg5Gd_4w' 0
    'cpd_Mg5Gd_E38413eV_cam3056x3056_pix2403nm_87R_2017_11003440_syn51_87R_Mg5Gd_4w' 0
    'cpd_Mg5Gd_E38413eV_cam3056x3056_pix2403nm_95R_2017_11003440_syn52_95R_Mg5Gd_8w' 0
    'cpd_Mg5Gd_E38402eV_cam5120x3840_pix1282nm_91R_2017_11003440_syn64_91R_Mg5Gd_4w' 0
    'cpd_Mg5Gd_E38402eV_cam5120x3840_pix1282nm_94R_2017_11003440_syn65_94R_Mg5Gd_8w' 0
    'cpd_Mg5Gd_E38390eV_cam5120x3840_pix1282nm_70L_2017_11003440_syn90_70L' 0
    'cpd_Mg5Gd_E38390eV_cam5120x3840_pix1282nm_72L_2017_11003440_syn91_72L' 0
    'cpd_Mg5Gd_E38390eV_cam5120x3840_pix1282nm_93R_2017_11003440_syn92_93R' 0
    'cpd_Mg5Gd_E38390eV_cam5120x3840_pix1282nm_76L_2017_11003440_syn93_76L' 0
    'cpd_Mg5Gd_E38390eV_cam5120x3840_pix1282nm_96R_2017_11003440_syn94_96R' 0
    'cpd_Mg5Gd_E38390eV_cam5120x3840_pix1282nm_69L_2017_11003440_syn95_69L' 0
    'cpd_Mg5Gd_E38390eV_cam5120x3840_pix1282nm_82L_2017_11003440_syn96_82L' 0
    'cpd_Mg5Gd_E38390eV_cam5120x3840_pix1282nm_82L_2017_11003440_syn96_82L_Mg5Gd_8w' 0
    'cpd_Mg5Gd_E38390eV_cam5120x3840_pix1282nm_74L_2017_11003440_syn97_74L' 0
    'cpd_Mg5Gd_E38390eV_cam5120x3840_pix1282nm_77L_2017_11003440_syn98_77L' 0
    'cpd_Mg5Gd_E28681eV_cam5120x3840_pix1280nm_62L_2017_11004016_syn001_62L_Mg5Gd_12w' 0
    'cpd_Mg5Gd_E33999eV_cam3056x3056_pix2420nm_96R_2018_11004263_syn004_96R_Mg5Gd_8w_a' 5
    'cpd_Mg5Gd_E33998eV_cam3056x3056_pix2420nm_96R_2018_11004263_syn004_96R_Mg5Gd_8w_b' 5
    'cpd_Mg5Gd_E33998eV_cam3056x3056_pix2420nm_76L_2018_11004263_syn007_76L_Mg5Gd_8w_a' 6
    'cpd_Mg5Gd_E33998eV_cam3056x3056_pix2420nm_76L_2018_11004263_syn007_76L_Mg5Gd_8w_b' 6
    'cpd_Mg5Gd_E33998eV_cam3056x3056_pix2420nm_93R_2018_11004263_syn008_93R_Mg5Gd_8w_a' 7
    'cpd_Mg5Gd_E33998eV_cam3056x3056_pix2420nm_93R_2018_11004263_syn008_93R_Mg5Gd_8w_b' 7
    'cpd_Mg5Gd_E33999eV_cam3056x3056_pix2420nm_64L_2018_11004263_syn009_64L_Mg5Gd_12w_a' 8
    'cpd_Mg5Gd_E33998eV_cam3056x3056_pix2420nm_64L_2018_11004263_syn009_64L_Mg5Gd_12w_b' 8
    'cpd_Mg5Gd_E33998eV_cam3056x3056_pix2420nm_90R_2018_11004263_syn011_90R_Mg5Gd_4w_a' 9
    'cpd_Mg5Gd_E33998eV_cam3056x3056_pix2420nm_90R_2018_11004263_syn011_90R_Mg5Gd_4w_b' 9
    'cpd_Mg5Gd_E33998eV_cam3056x3056_pix2420nm_85R_2018_11004263_syn014_85R_Mg5Gd_4w_a' 10
    'cpd_Mg5Gd_E33998eV_cam3056x3056_pix2420nm_85R_2018_11004263_syn014_85R_Mg5Gd_4w_b' 10
    'cpd_Mg5Gd_E33998eV_cam3056x3056_pix2420nm_97L_2018_11004263_syn015_97L_Mg5Gd_4w_a' 11
    'cpd_Mg5Gd_E33999eV_cam3056x3056_pix2420nm_97L_2018_11004263_syn015_97L_Mg5Gd_4w_b' 11
    'cpd_Mg5Gd_E33998eV_cam3056x3056_pix2420nm_91R_2018_11004263_syn022_91R_Mg5Gd_4w_a' 12
    'cpd_Mg5Gd_E33997eV_cam3056x3056_pix2420nm_91R_2018_11004263_syn022_91R_Mg5Gd_4w_b' 12
    'cpd_Mg5Gd_E33998eV_cam3056x3056_pix2420nm_72L_2018_11004263_syn023_72L_Mg5Gd_12w_a' 13
    'cpd_Mg5Gd_E33998eV_cam3056x3056_pix2420nm_72L_2018_11004263_syn023_72L_Mg5Gd_12w_b' 13
    'cpd_Mg5Gd_E33998eV_cam5120x3840_pix1281nm_72L_2018_11004263_syn032_72L_Mg5Gd_12w_kit' 0
    'cpd_Mg5Gd_E33998eV_cam5120x3840_pix1281nm_70L_2018_11004263_syn033_70L_Mg5Gd_12w_a' 14
    'cpd_Mg5Gd_E33999eV_cam5120x3840_pix1281nm_70L_2018_11004263_syn033_70L_Mg5Gd_12w_b' 14
    'cpd_Mg5Gd_E33998eV_cam5120x3840_pix1281nm_74L_2018_11004263_syn034_74L_Mg5Gd_8w_a' 15
    'cpd_Mg5Gd_E33999eV_cam5120x3840_pix1281nm_74L_2018_11004263_syn034_74L_Mg5Gd_8w_b' 15
    'cpd_Mg5Gd_E33998eV_cam5120x3840_pix1281nm_77L_2018_11004263_syn035_77L_Mg5Gd_8w_a' 16
    'cpd_Mg5Gd_E33999eV_cam5120x3840_pix1281nm_77L_2018_11004263_syn035_77L_Mg5Gd_8w_b' 16
    'cpd_Mg5Gd_E35000eV_cam3056x3056_pix2390nm_94R_2018_11005553_syn005_94R_Mg5Gd_8w' 0
    'cpd_Mg5Gd_E35000eV_cam5120x3840_pix1281nm_59R_2018_11005553_syn034_59R_Mg5Gd_12w' 0
    'cpd_Mg10Gd_E34014eV_cam3056x1621_pix2427nm_87L_2017_11003950_syn22_87L_Mg5Gd_4w_a' 17
    'cpd_Mg10Gd_E34014eV_cam3056x1621_pix2427nm_87L_2017_11003950_syn22_87L_Mg5Gd_4w_b' 17
    'cpd_Mg10Gd_E34276eV_cam3056x3056_pix2414nm_95L_2017_11003288_syn136_95L_Mg10Gd_8w' 0
    'cpd_Mg10Gd_E34276eV_cam3056x3056_pix2414nm_96L_2017_11003288_syn137_96L_Mg10Gd_8w' 0
    'cpd_Mg10Gd_E38374eV_cam3056x3056_pix2403nm_88L_2017_11003440_syn18_88L_Mg10Gd_4w' 0
    'cpd_Mg10Gd_E38374eV_cam3056x3056_pix2403nm_74R_2017_11003440_syn19_74R_Mg10Gd_8w' 0
    'cpd_Mg10Gd_E38374eV_cam3056x3056_pix2403nm_62R_2017_11003440_syn20_62R_Mg10Gd_12w' 0
    'cpd_Mg10Gd_E38374eV_cam3056x3056_pix2403nm_82R_2017_11003440_syn21_82R_Mg10Gd_8w' 0
    'cpd_Mg10Gd_E38402eV_cam3056x3056_pix2403nm_69R_2017_11003440_syn26_69R_Mg10Gd_12w' 0
    'cpd_Mg10Gd_E38402eV_cam3056x3056_pix2403nm_69R_2017_11003440_syn28_69R_Mg10Gd_12w' 0
    'cpd_Mg10Gd_E38402eV_cam3056x3056_pix2403nm_87L_2017_11003440_syn30_87L_Mg10Gd_4w' 0
    'cpd_Mg10Gd_E38402eV_cam3056x3056_pix2403nm_97R_2017_11003440_syn31_97R_Mg10Gd_4w' 0
    'cpd_Mg10Gd_E38402eV_cam3056x3056_pix2403nm_99R_2017_11003440_syn32_99R_Mg10Gd_4w' 0
    'cpd_Mg10Gd_E38402eV_cam3056x3056_pix2403nm_80R_2017_11003440_syn33_80R_Mg10Gd_8w' 0
    'cpd_Mg10Gd_E38402eV_cam3056x3056_pix2403nm_79R_2017_11003440_syn34_79R_Mg10Gd_8w' 0
    'cpd_Mg10Gd_E38402eV_cam3056x3056_pix2403nm_77R_2017_11003440_syn35_77R_Mg10Gd_8w' 0
    'cpd_Mg10Gd_E38402eV_cam3056x3056_pix2403nm_73R_2017_11003440_syn38_73R_Mg10Gd_8w' 0
    'cpd_Mg10Gd_E35000eV_cam5120x3840_pix1282nm_53L_2017_11003773_syn100_53L_Mg5Gd_4w' 0
    'cpd_Mg10Gd_E37000eV_cam3056x3056_pix2407nm_53L_2017_11003773_syn101_53L_Mg5Gd_4w' 0
    'cpd_Mg10Gd_E37000eV_cam3056x3056_pix2407nm_91L_2017_11003773_syn104_91L_Mg10Gd_4w' 0
    'cpd_Mg10Gd_E37000eV_cam3056x3056_pix2407nm_70R_2017_11003773_syn106_70R_Mg10Gd_12w' 0
    'cpd_Mg10Gd_E37000eV_cam3056x3056_pix2407nm_66R_2017_11003773_syn107_66R_Mg10Gd_12w' 0
    'cpd_Mg10Gd_E37000eV_cam3056x3056_pix2407nm_68R_2017_11003773_syn109_68R_Mg10Gd_12w' 0
    'cpd_Mg10Gd_E33998eV_cam3056x3056_pix2420nm_90L_2018_11004263_syn005_90L_Mg10Gd_4w_a' 18
    'cpd_Mg10Gd_E33998eV_cam3056x3056_pix2420nm_90L_2018_11004263_syn005_90L_Mg10Gd_4w_b' 18
    'cpd_Mg10Gd_E33998eV_cam3056x3056_pix2420nm_78R_2018_11004263_syn010_78R_Mg10Gd_8w_a' 19
    'cpd_Mg10Gd_E33998eV_cam3056x3056_pix2420nm_78R_2018_11004263_syn010_78R_Mg10Gd_8w_b' 19
    'cpd_Mg10Gd_E33998eV_cam3056x3056_pix2420nm_93L_2018_11004263_syn012_93L_Mg10Gd_8w_a' 20
    'cpd_Mg10Gd_E33998eV_cam3056x3056_pix2420nm_93L_2018_11004263_syn012_93L_Mg10Gd_8w_b' 20
    'cpd_Mg10Gd_E33998eV_cam3056x3056_pix2420nm_77R_2018_11004263_syn019_77R_Mg10Gd_8w_a' 21
    'cpd_Mg10Gd_E33998eV_cam3056x3056_pix2420nm_77R_2018_11004263_syn019_77R_Mg10Gd_8w_b' 21
    'cpd_Mg10Gd_E33998eV_cam3056x3056_pix2420nm_70R_2018_11004263_syn021_70R_Mg10Gd_12w_a' 22
    'cpd_Mg10Gd_E33998eV_cam3056x3056_pix2420nm_70R_2018_11004263_syn021_70R_Mg10Gd_12w_b' 22
    'cpd_Mg10Gd_E35000eV_cam3056x3056_pix2390nm_64R_2018_11005553_syn002_64R_Mg10Gd_12w' 0
    'cpd_Mg10Gd_E35000eV_cam5120x3840_pix1281nm_64R_2018_11005553_syn003_64R_Mg10Gd_12w_kit' 0
    'cpd_Mg10Gd_E35000eV_cam5120x3840_pix1281nm_60L_2018_11005553_syn028_60L_Mg10Gd_12w' 0
    'cpd_Mg10Gd_E35000eV_cam5120x3840_pix1281nm_68R_2018_11005553_syn033_68R_Mg10Gd_12w' 0
    'cpd_Mg10Gd_E35000eV_cam5120x3840_pix1281nm_56L_2018_11005553_syn035_56L_Mg10Gd_12w' 0
    'cpd_Ti_E45327eV_cam5120x3840_pix1280nm_39R_2019_11005842_syn027_18003_39R_bottom' -1
    'cpd_Ti_E45327eV_cam5120x3840_pix1280nm_39R_2019_11005842_syn027_18003_39R_middle' -1
    'cpd_Ti_E45327eV_cam5120x3840_pix1280nm_39R_2019_11005842_syn027_18003_39R_top' -1
    'cpd_Ti_E34005eV_cam3056x1201_pix2451nm_1R_2016_11001978_mah_35_1R_bottom' 23
    'cpd_Ti_E34005eV_cam3056x1201_pix2451nm_1R_2016_11001978_mah_36_1R_top' 23
    'cpd_Ti_E37000eV_cam3056x3056_pix2407nm_22L_2017_11003773_syn111_22L_Ti_8w' 0
    'cpd_Ti_E39486eV_cam5120x3840_pix1282nm_42L_2017_11003773_syn113_42L_Ti_12w' 0
    'cpd_Ti_E39486eV_cam5120x3840_pix1282nm_34R_2017_11003773_syn114_34R_Ti_8w' 0
    'cpd_Ti_E37000eV_cam3056x3056_pix2407nm_34R_2017_11003773_syn115_34R_Ti_8w' 0
    'cpd_Ti_E37000eV_cam3056x3056_pix2407nm_39R_2017_11003773_syn116_39R_Ti_12w' 0
    'cpd_Ti_E37000eV_cam3056x3056_pix2407nm_31R_2017_11003773_syn117_31R_Ti_8w' 0
    'cpd_Ti_E37000eV_cam3056x3056_pix2407nm_45L_2017_11003773_syn118_45L_Ti_12w' 0
    'cpd_Ti_E37000eV_cam3056x3056_pix2407nm_9L_2017_11003773_syn119_9L_Ti_4w' 0
    'cpd_Ti_E37000eV_cam3056x3056_pix2407nm_1L_2017_11003773_syn120_1L_Ti_4w' 0
    'cpd_Ti_E37000eV_cam3056x3056_pix2407nm_7L_2017_11003773_syn121_7L_Ti_4w' 0
    'cpd_Ti_E37000eV_cam3056x3056_pix2407nm_4R_2017_11003773_syn122_4R_Ti_4w' 0
    'cpd_Ti_E37000eV_cam3056x3056_pix2407nm_10L_2017_11003773_syn123_10L_Ti_4w' 0
    'cpd_Ti_E37000eV_cam3056x3056_pix2407nm_11L_2017_11003773_syn124_11L_Ti_4w' 0
    'cpd_Ti_E42000eV_cam5120x3840_pix1278nm_12L_2019_11005842_syn001_12l_ti_12w' 0
    'cpd_Ti_E45000eV_cam5120x3840_pix1279nm_12L_2019_11005842_syn003_12l_ti_12w' 0
    'cpd_Ti_E45000eV_cam5120x3840_pix1279nm_12L_2019_11005842_syn004_12l_ti_12w' 0
    'cpd_Ti_E45000eV_cam5120x3840_pix1280nm_41L_2019_11005842_syn005_41l_ti_12w' 0
    'cpd_Ti_E45000eV_cam5120x3840_pix1280nm_43L_2019_11005842_syn006_43l_ti_12w' 0
    'cpd_Ti_E45000eV_cam5120x3840_pix1280nm_46L_2019_11005842_syn007_46l_ti_12w' 0
    'cpd_Ti_E45000eV_cam5120x3840_pix1280nm_48R_2019_11005842_syn008_48r_ti_12w' 0 
    'cpd_Ti_E45000eV_cam5120x3840_pix1280nm_52R_2019_11005842_syn010_52r_ti_12w' 0
    'cpd_Ti_E45000eV_cam5120x3840_pix1280nm_25L_2019_11005842_syn011_25l_ti_8w' 0
    'cpd_Ti_E45000eV_cam5120x3840_pix1280nm_27L_2019_11005842_syn013_27l_ti_8w' 0
    'cpd_Ti_E45000eV_cam5120x3840_pix1280nm_28L_2019_11005842_syn014_28l_ti_8w' 0
    'cpd_Ti_E45000eV_cam5120x3840_pix1280nm_29L_2019_11005842_syn015_29l_ti_8w' 0
    'cpd_Ti_E45000eV_cam5120x3840_pix1280nm_30L_2019_11005842_syn016_30l_ti_8w' 0
    'cpd_Ti_E45000eV_cam5120x3840_pix1280nm_38R_2019_11005842_syn017_38r_ti_8w' 0
    'cpd_Ti_E45000eV_cam5120x3840_pix1280nm_11L_2019_11005842_syn018_11l_ti_4w' 0    
    'cpd_Ti_E45000eV_cam5120x3840_pix1280nm_14L_2019_11005842_syn019_14l_ti_4w' 0
    'cpd_Ti_E45000eV_cam5120x3840_pix1280nm_15L_2019_11005842_syn020_15l_ti_4w' 0
    'cpd_Ti_E45000eV_cam5120x3840_pix1280nm_16L_2019_11005842_syn021_16l_ti_4w' 0
    'cpd_Ti_E45000eV_cam5120x3840_pix1280nm_19L_2019_11005842_syn023_19l_ti_4w' 0
    'cpd_Ti_E45000eV_cam5120x3840_pix1280nm_19L_2019_11005842_syn025_19l_ti_4w_fastscan_stepwiserot' 0
    'cpd_Ti_E45000eV_cam5120x3840_pix1280nm_26L_2019_11005842_syn026_26l_ti_8w' 0 
    % long term samples
    'cpd_PEEK_E45000eV_cam5120x3840_pix1280nm_39L_2019_11005842_syn031_18005_39l_bottom' 0
    'cpd_PEEK_E45000eV_cam5120x3840_pix1280nm_39L_2019_11005842_syn031_18005_39l_middle' 0
    'cpd_PEEK_E45000eV_cam5120x3840_pix1280nm_39L_2019_11005842_syn031_18005_39l_top' 0
    'cpd_PEEK_E34005eV_cam3056x1701_pix2451nm_4L_2016_11001978_mah_20_4L_bottom' 25
    'cpd_PEEK_E34005eV_cam3056x1701_pix2451nm_4L_2016_11001978_mah_20_4L_top' 25
    'cpd_PEEK_E34005eV_cam3056x1501_pix2451nm_16R_2016_11001978_mah_27_16R_bottom' 26
    'cpd_PEEK_E34005eV_cam3056x1501_pix2451nm_16R_2016_11001978_mah_27_16R_top' 26
    'cpd_PEEK_E34005eV_cam3056x1501_pix2451nm_15R_2016_11001978_mah_28_15R_bottom'  27
    'cpd_PEEK_E34005eV_cam3056x1501_pix2451nm_15R_2016_11001978_mah_28_15R_top'  27
    'cpd_PEEK_E34005eV_cam3056x1201_pix2451nm_15R_2016_11001978_mah_32_15R_top_occd800_withoutpaper' -1
    'cpd_PEEK_E34005eV_cam3056x1201_pix2451nm_10R_2016_11001978_mah_37_10R_bottom' 28
    'cpd_PEEK_E34005eV_cam3056x1201_pix2451nm_10R_2016_11001978_mah_38_10R_top' 28
    'cpd_PEEK_E34005eV_cam3056x1201_pix2451nm_3L_2016_11001978_mah_39_3L_bottom' 29
    'cpd_PEEK_E34005eV_cam3056x1201_pix2451nm_3L_2016_11001978_mah_40_3L_top' 29
    'cpd_PEEK_E34005eV_cam3056x1201_pix2451nm_9R_2016_11001978_mah_41_9R_bottom' 30
    'cpd_PEEK_E34005eV_cam3056x1201_pix2451nm_9R_2016_11001978_mah_42_9R_top' 30
    'cpd_PEEK_E34000eV_cam3056x1801_pix2427nm_48L_2017_11003950_syn01_48L_PEEK_12w_b' 31
    'cpd_PEEK_E34000eV_cam3056x1801_pix2427nm_48L_2017_11003950_syn01_48L_PEEK_12w_c' 31
    'cpd_PEEK_E34000eV_cam3056x1801_pix2427nm_46R_2017_11003950_syn02_46R_PEEK_12w_a' 32
    'cpd_PEEK_E34000eV_cam3056x1801_pix2427nm_46R_2017_11003950_syn02_46R_PEEK_12w_b' 32
    'cpd_PEEK_E34000eV_cam3056x1801_pix2427nm_12R_2017_11003950_syn03_12R_PEEK_12w_a' 33
    'cpd_PEEK_E34001eV_cam3056x1801_pix2427nm_12R_2017_11003950_syn03_12R_PEEK_12w_b' 33
    'cpd_PEEK_E34001eV_cam3056x1801_pix2427nm_30R_2017_11003950_syn04_30R_PEEK_8w_a' 34
    'cpd_PEEK_E34001eV_cam3056x1801_pix2427nm_30R_2017_11003950_syn04_30R_PEEK_8w_b' 34
    'cpd_PEEK_E34001eV_cam3056x3056_pix2427nm_48L_2017_11003950_syn14_48L_PEEK_12w_a' 35
    'cpd_PEEK_E34001eV_cam3056x3056_pix2427nm_48L_2017_11003950_syn14_48L_PEEK_12w_b' 35
    'cpd_PEEK_E34001eV_cam3056x3056_pix2427nm_29R_2017_11003950_syn15_29R_PEEK_8w_a' 36
    'cpd_PEEK_E34001eV_cam3056x3056_pix2427nm_29R_2017_11003950_syn15_29R_PEEK_8w_b' 36
    'cpd_PEEK_E34001eV_cam3056x3056_pix2427nm_43R_2017_11003950_syn16_43R_PEEK_12w_a' 37
    'cpd_PEEK_E34001eV_cam3056x3056_pix2427nm_43R_2017_11003950_syn16_43R_PEEK_12w_b' 37
    'cpd_PEEK_E34014eV_cam3056x1621_pix2427nm_28R_2017_11003950_syn23_28R_PEEK_8w_a' 38
    'cpd_PEEK_E34014eV_cam3056x1621_pix2427nm_28R_2017_11003950_syn23_28R_PEEK_8w_b' 38
    'cpd_PEEK_E34276eV_cam3056x3056_pix2414nm_28R_2017_11003288_syn134_28R_PEEK_8w' 0
    'cpd_PEEK_E34276eV_cam3056x3056_pix2414nm_27R_2017_11003288_syn135_27R_PEEK_8w' 0
    'cpd_PEEK_E38290eV_cam5120x3840_pix1255nm_52L_2017_11003440_syn01_52L_PEEK_12w' 0
    'cpd_PEEK_E38290eV_cam3056x3056_pix2402nm_52L_2017_11003440_syn02_52L_PEEK_12w_ccd' 0
    'cpd_PEEK_E34166eV_cam5120x3840_pix1255nm_52L_2017_11003440_syn04_52L_PEEK_12w_dcm_cmos_recotest' -1
    'cpd_PEEK_E38234eV_cam5120x3840_pix1279nm_52L_2017_11003440_syn05_52L_PEEK_12w_dmm_cmos_recotest' -1
    'cpd_PEEK_E38234eV_cam3056x3056_pix2403nm_52L_2017_11003440_syn06_52L_PEEK_12w_dmm_ccd_recotest' -1
    'cpd_PEEK_E38374eV_cam3056x3056_pix2403nm_39L_2017_11003440_syn16_39L_PEEK_12w' 0
    'cpd_PEEK_E38374eV_cam3056x3056_pix2403nm_21R_2017_11003440_syn17_21R_PEEK_8w' 0
    'cpd_PEEK_E38402eV_cam3056x3056_pix2403nm_38L_2017_11003440_syn42_38L_PEEK_8w' 0
    'cpd_PEEK_E38402eV_cam3056x3056_pix2403nm_38L_2017_11003440_syn43_38L_PEEK_8w' 0
    'cpd_PEEK_E38390eV_cam5120x3840_pix1282nm_43R_2017_11003440_syn99_43R' 0
    'cpd_PEEK_E37000eV_cam3056x3056_pix2407nm_37L_2017_11003773_syn102_37L_PEEK_8w' 0
    'cpd_PEEK_E37000eV_cam3056x3056_pix2407nm_25R_2017_11003773_syn108_25R_PEEK_8w' 0
    'cpd_PEEK_E33998eV_cam3056x3056_pix2420nm_45R_2018_11004263_syn013_45R_PEEK_12w_a' 39
    'cpd_PEEK_E33999eV_cam3056x3056_pix2420nm_45R_2018_11004263_syn013_45R_PEEK_12w_b' 39
    'cpd_PEEK_E33999eV_cam3056x3056_pix2420nm_26R_2018_11004263_syn016_26R_PEEK_8w_a' 40
    'cpd_PEEK_E33998eV_cam3056x3056_pix2420nm_26R_2018_11004263_syn016_26R_PEEK_8w_b' 40
    'cpd_PEEK_E33999eV_cam3056x3056_pix2420nm_41R_2018_11004263_syn017_41R_PEEK_12w_a' 41
    'cpd_PEEK_E33998eV_cam3056x3056_pix2420nm_41R_2018_11004263_syn017_41R_PEEK_12w_b' 41
    'cpd_PEEK_E33998eV_cam3056x3056_pix2420nm_35L_2018_11004263_syn018_35L_PEEK_8w_a' 42
    'cpd_PEEK_E00000eV_cam3056x3056_pix2420nm_35L_2018_11004263_syn018_35L_PEEK_8w_c' 42
    'cpd_PEEK_E33998eV_cam3056x3056_pix2420nm_61L_2018_11004263_syn020_61L_PEEK_12w_a' 43
    'cpd_PEEK_E33998eV_cam3056x3056_pix2420nm_61L_2018_11004263_syn020_61L_PEEK_12w_b' 43
    'cpd_PEEK_E35000eV_cam5120x3840_pix1281nm_13R_2018_11005553_syn023_13R_PEEK_4w' 0
    'cpd_PEEK_E35000eV_cam5120x3840_pix1281nm_20R_2018_11005553_syn027_20R_PEEK_4w' 0
    'cpd_PEEK_E34005eV_cam3056x1501_pix2451nm_8L_2016_11001978_mah_26_8L_bottom' 44
    'cpd_PEEK_E34005eV_cam3056x1501_pix2451nm_8L_2016_11001978_mah_26_8L_top' 44
    'cpd_PEEK_E25000eV_cam3056x2201_pix1300nm_13R_2016_11001978_mah_10_13R_bottom' 45
    'cpd_PEEK_E25000eV_cam3056x2201_pix1300nm_13R_2016_11001978_mah_10_13R_top' 45
    'cpd_PEEK_E25000eV_cam3056x2201_pix1300nm_20R_2016_11001978_mah_11_20R_bottom' 24
    'cpd_PEEK_E25000eV_cam3056x2201_pix1300nm_20R_2016_11001978_mah_11_20R_top' 24
    'cpd_Mg10Gd_E38402eV_cam3056x3056_pix2403nm_100AL_2017_11003440_syn47_100AL_Mg5Gd_4w' 0
    };
fprintf( '\n number of cpd scans : %u', size( cpd_list, 1 ) )
%% Check if new recos arrived
cprintf( 'Red', '\nTest for new scans:' )
s = scans.cpd.all;
mask = zeros( [numel( s ),1], 'logical' );
for n = 1:numel( s )
    sn = s(n);
    for m = 1:size( cpd_list, 1 )
        cpdm = cpd_list{m,1};
        if strcmp( sn.name, cpdm )
            mask(n) = 1;
        end
    end
end
s(mask) = [];
if ~isempty( s )
    fprintf( '\n' )
    fprintf( '\n''%s''', s(:).name )
    fprintf( '\n\n' )
    error( 'ADD SCANS TO LIST' )
else
    fprintf( ' NONE' )
end

%% Test for obsolete scans
cprintf( 'Red', '\nTest for obsolete scans:' )
s = scans.cpd.all;
c = cpd_list;
mask = zeros( [size( c, 1 ),1], 'logical' );
for n = 1:size( c, 1 )
    cname = c{n,1};
    for m = 1:numel( s )
        sname = s(m).name;
        if strcmp( sname, cname )
            mask(n) = 1;
        end
    end
end
o = c(~mask);
if ~isempty( o )
    fprintf( '\n' )
    fprintf( '\n''%s''', o{:} )
    fprintf( '\n\n' )
    error( 'REMOVE SCANS FROM LIST' )
else
    fprintf( ' NONE' )
end

%% Plot pixel size binned
s = scans.cpd.all;
name = 'Effective pixel size binned';
f0 = figure('Name', name );
plot([s.effective_pixel_size_binned],'.')
title( name )
xlabel( 'scan number' )
ylabel( 'effective pixelsize binned' )
axis tight
drawnow
fig_filename = sprintf( '%s%s.png', mfile_path, regexprep( f0.Name, '\ |:', '_') ) ;
saveas( f0, fig_filename );
%% LIST : Stitch, Single, Other
scan_type = [cpd_list{:,2}];
mask = scan_type > 0;
scans_stitch = cpd_list(mask,1);
fprintf( '\n number of cpd scans to stitch : %u', numel( scans_stitch ) )
mask = scan_type == 0;
scans_single_height = cpd_list(mask,1);
fprintf( '\n number of cpd single height : %u', numel( scans_single_height ) )
mask = scan_type == -1;
scans_other = cpd_list(mask,1);
fprintf( '\n number of cpd other : %u', numel( scans_other ) )

%% Stitch struct
cprintf( 'red', '\n\nCreate stitch struct' )
scans_stitch_struct(max( scan_type )).scans = {};
for nn = 1:max( scan_type )
    mask = scan_type == nn;
    scans_stitch_struct(nn).scans = cpd_list( mask )';
    scans_stitch_struct(nn).num_scans = sum( mask(:) );
    for mm = 1:scans_stitch_struct(nn).num_scans
        name = scans_stitch_struct(nn).scans{mm};
        match = strcmp( name, {scans_proc.name} );
        scans_stitch_struct(nn).scan_structs(mm) = scans_proc(match);
        
        % Create string for height level
        str = scans_stitch_struct(nn).scan_structs(mm).folder_name;
        c = regexp(str, '_', 'split');
        % Delete height level indicator
        c(end) = [];
        if mm == 1
            c10 = c{1};
            c20 = c{2};
        else
            c11 = c{1};
            c21 = c{2};
            if strcmp( c10, c11 )
                cs1 = c10;
            else
                cs1 = [ c10 '_' c11 ];
            end
            if strcmp( c20, c21 )
                cs2 = c20;
            else
                cs2 = [ c20 '_' c21 ];
            end
        end
    end
    folder_name = [ sprintf( '%s_%s', cs1, cs2 ) sprintf( '_%s', c{3:end} ) ];
    s = scans_stitch_struct(nn).scan_structs(1);
    name = sprintf( '%s_%s_E%05.0feV_cam%ux%u_pix%04.0fnm%s_%u_%u_%s', ...
        s.type, ...
        s.material, ...
        s.energy, ...
        s.im_shape_raw, ...
        1000*s.effective_pixel_size, ...
        s.sample, ...
        s.beamtime_year, ...
        s.beamtime_id, ...
        folder_name);
    scans_stitch_struct(nn).name = name;
    fprintf( '\n %s', name )
    if scans_stitch_struct(nn).num_scans ~= numel( scans_stitch_struct(nn).scan_structs )
        error( '\nCheck scans structs\n' )
    end
end
fprintf( '\n #elements: %u', numel( scans_stitch_struct ) )

%% Single height struct
cprintf( 'red', '\n\nCreate single height struct' )
for nn = numel(scans_single_height):-1:1
    name = scans_single_height{nn};
    match = strcmp( name, {scans_proc.name} );
    single_height_scans_struct(nn) = scans_proc(match);
end
fprintf( '\n #elements: %u', numel( single_height_scans_struct ) )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% RESAMPLIG %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cprintf( 'red', '\n\nSTART RESAMPLING PROCEDURE' )
%% Height scans %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
total_number_scans = numel( scans_stitch_struct ) + numel( single_height_scans_struct );
if ~exist( 'total_scan_counter', 'var' )
    total_scan_counter = 0;
end
if run_height_scans
    cprintf( 'red', '\n Height scans to stitch' )
    for nn = 1:numel( scans_stitch_struct )
        scan = scans_stitch_struct(nn);
        name = scan.name;
        fprintf( '\n\n%2u: %s', nn, name )
        fprintf( '\n height levels: %u', scan.num_scans )
        
        try
            stitch_struct = stitch_volumes_2( scan );
        catch ME
            keyboard
        end
        
        effective_pixel_size_binned = stitch_struct(1).effective_pixel_size_binned;
        
        % Resample
        fprintf( '\n Resampling.' )
        
        % Concatenate volume
        V = cat(3, stitch_struct(stitch_struct(1).stitch_order).vol );
        zval = cat( 2, stitch_struct(stitch_struct(1).stitch_order).zval_stitch);
        clear stitch_struct
        
        % Original grid
        [x,y,~] = size( V );
        xx = 1:x;
        yy = 1:y;
        zz = zval / effective_pixel_size_binned;
        [X,Y,Z] = meshgrid( xx, yy, zz );
        fprintf( '\n Sample point ranges:' )
        fprintf( '\n  x: %f %f', xx(1), xx(end) )
        fprintf( '\n  x: %f %f', yy(1), yy(end) )
        fprintf( '\n  x: %f %f', zz(1), zz(end) )
        
        % Query grid
        k = 5 / effective_pixel_size_binned;
        k0 = max( k, 1 );
        xq = k0:k:x;
        yq = k0:k:y;
        zq = zz(1):k:zz(end);
        [Xq,Yq,Zq] = meshgrid( xq, yq, zq );
        fprintf( '\n Query point ranges:' )
        fprintf( '\n  x: %f %f', xq(1), xq(end) )
        fprintf( '\n  y: %f %f', yq(1), yq(end) )
        fprintf( '\n  z: %f %f', zq(1), zq(end) )
        
        %Vq = interp3( vol, Xq, Yq, Zq );
        Vq = interp3( X, Y, Z, V, Xq, Yq, Zq );
        clear X Y Z V Xq Yq Zq;
        [x,y,z] = size( Vq );
        
        % Check volume
        if sum( isinf( Vq(:) ) ) || sum( isnan( Vq(:) ) )
            error( 'Check volume' )
        end
        vmin = min( Vq(:) );
        vmax = max( Vq(:) );
        vmean = mean( Vq(:) );
        vstd = std( Vq(:) );
        scan.vmin = vmin;
        scan.vmax = vmax;
        scan.vmean = vmean;
        scan.vstd = vstd;
        fprintf( '\n [min, max, max-min, mean, std] = [%f, %f, %f, %f, %f] ', vmin, vmax, vmax - vmin, vmean, vstd )
        %[N,edges] = histcounts( Vq(:), 1024 );
        [N,edges] = histcounts( Vq(:), 1024, 'BinLimits', bin_limits );
        scan.histcounts.N = N;
        scan.histcounts.edges = edges;
        total_scan_counter = total_scan_counter + 1;
        scans_resampled( total_scan_counter ) = scan;
        scans_resampled( total_scan_counter ).name = name;
        
        %% Save resampled volume
        fprintf( '\n Saving.' )
        save_path = sprintf( '%s%s', parpath, name );
        CheckAndMakePath( save_path )
        parfor kk = 1:size( Vq, 3 )
            filename = sprintf( '%s/vol_%06u.tif', save_path, kk );
            write32bitTIFfromSingle( filename, Vq(:,:,kk) );
        end
        
        %% Save ortho slices
        xx = round( x / 2 );
        yy = round( y / 2 );
        zz = round( z / 2 );
        % Save ortho slices x
        im = rot90( FilterOutlier( squeeze( Vq(xx,:,:) ), 0.02, '', 0, 0 ), -1 );
        filename = sprintf( '%s/%s_x%06uof%06u.png', overviewpath, name, xx, x );
        imwrite( normat(im), filename );
        % Save ortho slices y
        im = rot90( FilterOutlier( squeeze( Vq(:,yy,:) ), 0.02, '', 0, 0 ),-1 );
        filename = sprintf( '%s/%s_y%06uof%06u.png', overviewpath, name, yy, y );
        imwrite( normat(im), filename );
        % Save ortho slices z
        im = FilterOutlier( squeeze( Vq(:,:,zz) ), 0.02, '', 0, 0 );
        filename = sprintf( '%s/%s_z%06uof%06u.png', overviewpath, name, zz, z );
        imwrite( normat(im), filename );
        clear Vq
    end
    filename = [ mfile_path 'scans_resampled.mat'];
    save( filename, 'scans_resampled' );
else
    total_scan_counter = numel( scans_stitch_struct );
    load('scans_resampled.mat');
end
%% Single height scans %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ca
if run_single_height_scans
    cprintf( 'red', '\n\nSINGLE HEIGHT SCANS' )
    single_height_scans_to_process = [];
    if isempty( single_height_scans_to_process )
        single_height_scans_to_process = single_height_scans_struct;
    end
    %single_height_scans_to_process = single_height_scans_struct(end);
    
    for nn = 17:numel( single_height_scans_to_process )
        scan = single_height_scans_to_process(nn);
        name = scan.name;
        
        fprintf( '\n%s', name )
        data_path = scan.data_path;
        bin = double( scan.bin );
        effective_pixel_size = scan.effective_pixel_size;
        effective_pixel_size_binned = bin * effective_pixel_size;
        
        % Print parameters
        fprintf( '\n  bin: %u', bin )
        fprintf( '\n  effective pixelsizse: %f', effective_pixel_size )
        fprintf( '\n  effective pixelsizse_binned: %f', effective_pixel_size_binned )
        
        % Read volume
        im_struct = dir( [ data_path filesep '*.tif' ] );
        vol = imread( [data_path filesep im_struct(1).name] );
        vol = zeros( [size(vol) numel( im_struct )], 'single' );
        fprintf( '\n Reading.' )
        parfor kk = 1:numel( im_struct )
            filename = [data_path filesep im_struct(kk).name ];
            vol(:,:,kk) = imread( filename, 'tif' );
        end
        [x,y,z] = size( vol );
        fprintf( ' size: %u %u %u', x, y, z )
        
        % Resample
        fprintf( '\n Meshgrid' )
        k = 5 / effective_pixel_size_binned;
        k0 = max( k, 1 );
        xq = k0:k:x;
        yq = k0:k:y;
        zq = k0:k:z;
        [Xq,Yq,Zq] = meshgrid( xq, yq, zq );
        fprintf( '\n Query point ranges:' )
        fprintf( '\n  x: %f %f', xq(1), xq(end) )
        fprintf( '\n  y: %f %f', yq(1), yq(end) )
        fprintf( '\n  z: %f %f', zq(1), zq(end) )
        fprintf( '\n Interpolate.' )
        Vq = interp3( vol, Xq, Yq, Zq );
        [x,y,z] = size( Vq );
        fprintf( ' resampled size: %u %u %u', x, y, z )
        
        % Check volume
        if sum( isinf( Vq(:) ) ) || sum( isnan( Vq(:) ) )
            error( 'Check volume' )
        end
        vmin = min( Vq(:) );
        vmax = max( Vq(:) );
        vmean = mean( Vq(:) );
        vstd = std( Vq(:) );
        total_scan_counter = total_scan_counter + 1;
        scans_resampled(total_scan_counter).scans = scan;
        scans_resampled(total_scan_counter).name = name;
        scans_resampled(total_scan_counter).vmin = vmin;
        scans_resampled(total_scan_counter).vmax = vmax;
        scans_resampled(total_scan_counter).vmean = vmean;
        scans_resampled(total_scan_counter).vstd = vstd;
        fprintf( '\n [min, max, max-min, mean, std] = [%f, %f, %f, %f, %f] ', vmin, vmax, vmax - vmin, vmean, vstd )
        [N,edges] = histcounts( Vq(:), 1024, 'BinLimits', bin_limits );
        scans_resampled(total_scan_counter).histcounts.N = N;
        scans_resampled(total_scan_counter).histcounts.edges = edges;
        
        
        %% Save resampled volume
        fprintf( '\n Saving.' )
        save_path = sprintf( '%s%s', parpath, name );
        CheckAndMakePath( save_path )
        fprintf( '\n save_path: %s', save_path );
        parfor kk = 1:size( Vq, 3 )
            filename = sprintf( '%s/vol_%06u.tif', save_path, kk );
            write32bitTIFfromSingle( filename, Vq(:,:,kk) );
        end
        
        %% Save ortho slices
        xx = round( x / 2 );
        yy = round( y / 2 );
        zz = round( z / 2 );
        % Save ortho slices x
        im = rot90( FilterOutlier( squeeze( Vq(xx,:,:) ), 0.02, '', 0, 0 ), -1 );
        filename = sprintf( '%s/%s_x%06uof%06u.png', overviewpath, name, xx, x );
        imwrite( normat(im), filename );
        % Save ortho slices y
        im = rot90( FilterOutlier( squeeze( Vq(:,yy,:) ), 0.02, '', 0, 0 ),-1 );
        filename = sprintf( '%s/%s_y%06uof%06u.png', overviewpath, name, yy, y );
        imwrite( normat(im), filename );
        % Save ortho slices z
        im = FilterOutlier( squeeze( Vq(:,:,zz) ), 0.02, '', 0, 0 );
        filename = sprintf( '%s/%s_z%06uof%06u.png', overviewpath, name, zz, z );
        imwrite( normat(im), filename );
    end
end

if total_scan_counter ~= total_number_scans
    error( '\n Something is fishy. Check and optionally save ''scans_resampled''.' )
else
    filename = [ mfile_path 'scans_resampled.mat'];
    save( filename, 'scans_resampled' );
end
%% Figures %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ind = 1:numel( scans_resampled );
name = 'Volume: min max mean std';
if ~ishandle( f1 )
    f1 = figure( 'Name', name );            
end
vmin = [scans_resampled(ind).vmin];
vmax = [scans_resampled(ind).vmax];
vmean = [scans_resampled(ind).vmean];
vstd = [scans_resampled(ind).vstd];
Y = [ vmin' vmax' vmean' vstd'];
plot( Y , '.')
title( name )
xlabel( 'scan number' )
legend( { 'min', 'max', 'mean', 'std' } )
axis tight
drawnow
fig_filename = sprintf( '%s%s.png', mfile_path, regexprep( name, '\ |:', '_') ) ;
saveas( f1, fig_filename );
%% Figure: histogram
name = 'Histograms';
if ~exist( 'f2', 'var' ) || ~ishandle( f2 )
    f2 = figure( 'Name', name );            
end
num_scans = numel( scans_resampled );
num_bins = numel( scans_resampled(1).histcounts.N );
X = meshgrid( 1:num_scans, 1:num_bins );
clear Y Z
for nn = num_scans:-1:1
    b = scans_resampled(nn).histcounts.edges;
    Y(:,nn) = ( b(2:end) + b(1:end-1) ) / 2;
    Z(:,nn) = scans_resampled(nn).histcounts.N;
end
whos X Y Z
plot3( X, Y, Z )
title( name )
xlabel( 'scan number' )
ylabel( 'bins' )
zlabel( 'counts' )
axis tight
drawnow
fig_filename = sprintf( '%s%s.png', mfile_path, regexprep( name, '\ |:', '_') ) ;
saveas( f2, fig_filename );
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf( '\nFINISHED.\n' )
