%ca
%clear all

read_scans = 0;

%% Beamtimes
beamtime = { ...
    2016 11001978
    2017 11003950
    2017 11004016
    2017 11003288
    2017 11003440
    2017 11003440
    2018 11004263
    2018 11004936
    2018 11005553};

%% Output directory
parpath = '/asap3/petra3/gpfs/p05/2018/data/11004936/processed/segmentation/';
overviewpath = [parpath 'overview'];

%% Read scans %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if read_scans
    
    scans =[];
    n_cpd = 0;
    n_load = 0;
    n_other = 0;
    
    % Loop over beamtimes
    for nn = 1:length( beamtime )
        
        % beamtime folder
        beamtime_year = beamtime{nn,1};
        beamtime_id = beamtime{nn,2};
        proc = sprintf( '/asap3/petra3/gpfs/p05/%u/data/%u/processed', beamtime_year, beamtime_id );
        
        % scan folder
        scan_struct = dir( [proc '/*Mg*'] );
        
        fprintf( '\nexp %2u: year %u, beamtime %u, %u folders with Mg', nn, beamtime_year, beamtime_id, length( scan_struct ) );
        
        % Loop over scans
        for mm = 1:length( scan_struct )
            scan_name = scan_struct(mm).name;
            scan_path = [ scan_struct(mm).folder filesep scan_name ];
            
            if isdir( scan_path )
                
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
                    vol_name = floats_namecell{1};
                    vol_path = [reco_path filesep vol_name];
                    
                    %% Read log file
                    log_file = [reco_path '/reco.log'];
                    fid = fopen( log_file );
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
                    
                    %% Binning factor
                    t = regexp( vol_name, 'bin|Bin' );
                    if ~isempty( t )
                        bin = str2double( vol_name(t+3) );
                    else
                        bin = 1;
                    end
                    pixel_size = effective_pixel_size * bin;
                    
                    
                    %% Unique name
                    unique_scan_name = sprintf( '%u_%u_%s', beamtime_year, beamtime_id, scan_name);
                    
                    %% Create parameter struct
                    tmp.unique_name =unique_scan_name;
                    tmp.bin = bin;
                    tmp.effective_pixel_size = effective_pixel_size;
                    tmp.full_reco_path = vol_path;
                    
                    %% Sort scan type
                    expr = '_\d\d\d|load|test|push|drill';
                    t = regexp( scan_name, expr, 'once');
                    
                    if isempty( t )
                        % CPD
                        n_cpd = n_cpd + 1;
                        scans.cpd(n_cpd) = tmp;
                        tag = 'c';
                    else
                        
                        expr = 'test' ;
                        t = regexp( scan_name, expr );
                        if isempty( t )
                            % load
                            n_load = n_load + 1;
                            scans.load(n_load) = tmp;
                            tag = 'l';
                        else
                            % rest
                            n_other = n_other + 1;
                            scans.other(n_other) = tmp;
                            tag = 'o';
                        end
                    end
                    
                    fprintf( '\n %2u. bin: %u, pixelsize: %f, tag: %s, name: %s', mm, bin, effective_pixel_size, tag, unique_scan_name )
                    
                end
            end
        end
        
        
        
    end
    
    %% Save parameter struct
    filename = [ parpath 'scan_paramter_struct.m'];
    save( filename, 'scans' );
    
    fprintf( '\nRead scans finished.\n' )
end

%% Stitch and resample %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
process.tags = {'c'};
resample_and_stitch = 1;
if resample_and_stitch
    fprintf( '\n\nResampling and stitching' )
    
    %% Loop over scan types
    for nn = 1:numel( process.tags )
        
        tag = process.tags{nn};
        
        fprintf( '\n Scan type: %s', tag );
        
        switch tag
            case 'c'
                ss = scans.cpd;
                outpath = [ parpath 'cpd/' ];
            case 'l'
                ss = scans.load;
                outpath = [ parpath 'load/' ];
            case 'o'
                ss = scans.other;
                outpath = [ parpath 'other/' ];
        end
        
        CheckAndMakePath( parpath )
        CheckAndMakePath( overviewpath );
        
        %% Loop over scans
        for mm = 1:numel( ss )
            
            s = ss(mm);
            
            fprintf( '\n %s', s.unique_name )
            
            
            
            
        end
    end
end
% %% Scan folders
% filename = [ parpath 'foldernames.txt'];
% fid = fopen( filename, 'w' );
% for ss = 1:numel( scans )
%     fprintf( fid, '%s\n', scans{ss} );
% end
% fclose( fid );


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf( '\nFINISHED.\n' )
