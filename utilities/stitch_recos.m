function stitch_recos( scan_path )
% Stich reconstructed volumes using log file information.
%
%
% Written by J. Moosmann 

if nargin < 1 
    scan_path = '/asap3/petra3/gpfs/p05/2018/data/11004263/processed/syn004_96R_Mg5Gd_8w';
end
if nargin < 2
    reco_subfolder = 'float_rawBin2';
end

%% Main %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Scans to stitch
scan_struct = dir( [scan_path '*' ] );

num_scans = numel( scan_struct );

% loop over scans
for nn = 1:num_scans
    
    % Name and folder
    name = scan_struct(nn).name;
    full_path = [scan_struct(nn).folder filesep name];
    
    %% Reco log
    reco_log = [full_path '/reco/reco.log' ];
    if ~exist( reco_log, 'file' )
        fprintf( '\nReco log not found!\n' )
        break
    end
    fid = fopen( reco_log );
    c = textscan( fid, '%s', 'Delimiter', {'\n', '\r'} );
    c = c{1};
    fclose( fid );
    
    %% Projection shape
    for ll = 1:numel( c )
        t = regexp( c{ll}, 'im_shape_raw' );
        if t
            cc = textscan( c{ll}, '%*s : %u %u');
            im_shape_raw = [ cc{1} cc{2} ];
            break
        end
    end
    
    %% Binning factor    
    for ll = 1:numel( c )
        t = regexp( c{ll}, 'raw_binning_factor' );
        if t
            cc = textscan( c{ll}, '%*s : %u %u');
            bin = [ cc{1} cc{2} ];
            break
        end
    end
        
    %% ROI
    for ll = 1:numel( c )
        t = regexp( c{ll}, 'raw_roi' );
        if t
            cc = textscan( c{ll}, '%*s : %u %u');
            raw_roi = [ cc{1} cc{2} ];
            break
        end
    end
    
    %% Effective pixelsize
    for ll = 1:numel( c )
        t = regexp( c{ll}, 'effective_pixel_size' );
        if t
            cc = textscan( c{ll}, '%*s : %f %*s');
            effective_pixel_size = cc{1};
            break
        end
    end
    
    %% Scan log
    scan_log = [regexprep( full_path, 'processed', 'raw' ) filesep name 'scan.log'];
     if ~exist( scan_log, 'file' )
        fprintf( '\nScan log not found!\n' )
        break
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
            z = str2double( cc{1}{2} );
            break
        end
    end
    
    fprintf( '\n %s, raw_shape: %u %u, raw_roi: %u %u, bin: %u, z: %f, pixelsizse: %f', name, im_shape_raw, raw_roi, bin, z, effective_pixel_size )
    
    
    %% Read volume
    vol_path = [full_path '/reco/' reco_subfolder];
    im_struct = dir( [ vol_path filesep '*.tif' ] );
    vol = imread( [vol_path filesep im_struct(1).name] );
    vol = zeros( [size(vol) numel( im_struct ) num_scans], 'single' );
    fprintf( ' Reading.' )
%     parfor kk = 1:numel( im_struct )
%         filename = [vol_path filesep im_struct(kk).name ];
%         vol(:,:,nn) = imread( filename, 'tif' );
%     end
%     
    
        
    
end




fprintf( '\n' )
