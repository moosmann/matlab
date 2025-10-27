function [s, vol] = stitch_volumes( scan_path, scan_subfolder, reco_subfolder, stitched_volume_path, scan_mask, noisecut, testing )
% Stich reconstructed volumes using log file information.
%
% ARGUMENTS
% scan_path : string. absolute path PATTERN common for all height levels or
%   cell of string. Note that default output paths differ if cell of
%   strings is used.
% scan_subfolder : string. default: 'reco'. typically 'reco' or 'reco_phase'
% reco_subfolder : string. default: 'float_rawBin2'. Subfolder to 'reco' to
%   be used for stitching.
% stitched_volume_path : string. default: scan_path (without trailing
%   underscores) or, if scan_path is a cell, the last cell string + 'stitched_volume'
% scan_mask : binary mask. Choose a subset of the scans found to be
% stitchted using a binary mask with a length of the number of scans found.
% noisecut : string = 'new' (default) or 'old'. old methods works for certain scans
%   only. new methods checks SNR by slicewise calculation of std/mean and
%   matches the SNR of the two volumes within the overlap region.
% testing: bool. default 0. if 1 data is neihter read or written. Used for
%   testing if scripts runs at all.
%
% RETRURNS
% s : struct containing individual arrays and full information for
%   stitching
% vol : stitched volume array
%
% Written by J. Moosmann
dbstop if error

%% ARGUMENTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
if nargin < 1
    scan_path  = ...
        {'/asap3/petra3/gpfs/p07/2025/data/11022007/processed/00015_cdma_lnkn1_a',...
        '/asap3/petra3/gpfs/p07/2025/data/11022007/processed/00016_cdma_lnkn1_b',...
        '/asap3/petra3/gpfs/p07/2025/data/11022007/processed/00017_cdma_lnkn1_c'};
   % scan_path =
   % '/asap3/petra3/gpfs/p07/2023/data/11017206/processed/itaw012_cet548a_OO01_Oo';
end
if nargin < 2
    %scan_subfolder = 'reco_phase/tie_regPar1p00';
    scan_subfolder = 'reco';
end
if nargin < 3
    reco_subfolder = 'reco2/float_rawBin2';
end
if nargin < 4
    stitched_volume_path = ...
    '/asap3/petra3/gpfs/p07/2025/data/11022007/processed/00015_cdma_lnkn1';
end
if nargin < 5
    scan_mask = [];
    %scan_mask = [1 1 0 0 0];
end
if nargin < 6
    noisecut = 'none';
    %noisecut = 'new';
    %noisecut = 'old';
end
if nargin < 7
    testing = 0;
end
% stitch level factor for old noisecut procedure
stitch_level_fac = 1.5;
normalize_overlap = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Output path
if isempty(stitched_volume_path)
    % Remove trailing underscores
    if iscell(stitched_volume_path)
        stitched_volume_path = scan_path{1};
    else
        stitched_volume_path = scan_path;
    end
    while strcmp(stitched_volume_path(end), '_')
        stitched_volume_path(end) = [];
    end
    [~,stitched_name]=fileparts(stitched_volume_path);
    stitched_volume_path  = sprintf( '%s/%s/%s', stitched_volume_path, scan_subfolder, reco_subfolder);
    stitched_volume_figure  = sprintf( '%s/%s/stitched_volume_%s', stitched_volume_path, scan_subfolder, reco_subfolder);
    stitched_volume_log_path  = sprintf( '%s/%s', stitched_volume_path, scan_subfolder);    
else
    [~,stitched_name]=fileparts(stitched_volume_path);
    CheckAndMakePath( stitched_volume_path )
    stitched_volume_path = [stitched_volume_path filesep scan_subfolder];
    CheckAndMakePath( stitched_volume_path )
    stitched_volume_figure  = sprintf( '%s/stitched_volume', stitched_volume_path);
    stitched_volume_log_path  = stitched_volume_path;
    stitched_volume_path = [stitched_volume_path filesep reco_subfolder];   
end
CheckAndMakePath( stitched_volume_path )
fprintf('\noutpath: %s' , stitched_volume_path)
fprintf('\nstitched name: %s' , stitched_name)
ca;
OpenParpool(0.5,0,stitched_volume_path,0,40);
%% Read parameters and data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('\nRead parameters and volumes')
if ~iscell( scan_path )
    fprintf( '\nscan_path: %s', scan_path )
    % Scans to stitch
    scan_struct = dir( [scan_path '*' ] );
    [~, scan_name ] = fileparts( scan_path );
    mask = ~strcmp( {scan_struct.name}, scan_name );
    scan_struct = scan_struct( mask );
    if ~isempty(scan_mask)
        scan_struct = scan_struct( logical( scan_mask ) );
    end
    % Allocate containing volumes
else
    for nn = numel( scan_path ):-1:1
        scan_struct(nn).full_path = scan_path{nn};
        [folder, name] = fileparts( scan_path{nn});
        scan_struct(nn).name = name;
        scan_struct(nn).folder = folder;
    end
end
num_scans = numel( scan_struct );

% loop over scans
s(num_scans).full_path = [];
fprintf('\nScans found: %u',num_scans)
fprintf('\nScans to be stitched:')
for nn = 1:num_scans
    fprintf('\n%s', scan_struct(nn).name)
end
fprintf('\nStart stitching:')
for nn = 1:num_scans
    
    % Name and folder
    name = scan_struct(nn).name;
    full_path = [scan_struct(nn).folder filesep name];
    s(nn).full_path = full_path;
    
    % Reco log
    reco_log = [full_path filesep scan_subfolder filesep 'reco*.log' ];
    rls = dir( reco_log );
    if ~isempty( rls )
        reco_log = [rls(end).folder filesep rls(end).name];
    else
        fprintf( '\nReco log not found!\n' )
        break
    end
    copy_log_fn = sprintf( '%s/%s_%s', stitched_volume_log_path, name, rls(end).name);
    if ~exist( copy_log_fn, 'file')
        copyfile( reco_log, copy_log_fn )
    end
    fid = fopen( reco_log );
    c = textscan( fid, '%s', 'Delimiter', {'\n', '\r'} );
    c = c{1};
    fclose( fid );
    s(nn).reco_log = reco_log;
    
    % Projection shape
    for ll = 1:numel( c )
        %t = regexp( c{ll}, 'im_shape_raw' );
        t = regexp( c{ll}, {'im_shape_raw', 'raw_image_shape', 'raw_image_shape_raw'} );
        %if t
        if sum( [t{:}] )
            cc = textscan( c{ll}, '%*s : %u %u');
            %im_shape_raw = [ cc{1} cc{2} ];
            if  isempty(cc{1})
                cc = textscan( c{ll}, '%*s [%u %u]');
            end
            im_shape_raw = [cc{:}];
            break
        end
    end
    im_shape_raw = double(im_shape_raw);
    s(nn).im_shape_raw = im_shape_raw;
    
    % Binning factor
    for ll = 1:numel( c )
        % t = regexp( c{ll}, 'raw_binning_factor' );
        t = regexp( c{ll}, {'raw_bin', 'raw_binning_factor', 'raw_binning'} );
        %if t
        if sum( [t{:}] )
            cc = textscan( c{ll}, '%*s : %u');            
            bin = cc{:};
            if isempty(bin)
                cc = textscan( c{ll}, '%*s %u');
                bin = cc{:};
            end
            if ~isempty(bin)
                break
            end
        end
    end
    bin = double(bin);
    s(nn).bin = bin;
    
    % ROI
    for ll = 1:numel( c )
        t = regexp( c{ll}, 'raw_roi' );
        if t
            cc = textscan( c{ll}, '%*s : %u %u');
            raw_roi = [ cc{1} cc{2} ];
            if isempty(cc)
                cc = textscan(c{ll},'%*s %u %u');
            end
            if isempty(cc)
                cc = textscan(c{ll},'%*s %s');
                isequal(cc{1}{1},'[]');
                raw_roi = [];
            end
            break
        end
    end
    if ~exist( 'raw_roi', 'var' ) || ( exist( 'raw_roi', 'var' ) && isempty( raw_roi ) )
        raw_roi(1) = 1;
        raw_roi(2) = im_shape_raw(2);
    end
    s(nn).raw_roi = raw_roi;
    
    % Effective pixelsize
    for ll = 1:numel( c )
        t = regexp( c{ll}, {'effective_pixel_size', 'eff_pixel_size', 'par.eff_pixel_size'} );
        if sum( [t{:}] )
            cc = textscan( c{ll}, '%*s : %f %*s');
            effective_pixel_size = cc{1};
            break
        end
    end
    s(nn).effective_pixel_size = effective_pixel_size;
    s(nn).effective_pixel_size_binned = effective_pixel_size * double( bin );

    % Scan log
    s_stage_z = [];
    scan_log = [regexprep( full_path, 'processed', 'raw' ) filesep name 'scan.log'];
    if ~exist( scan_log, 'file' )
        d = dir([regexprep( full_path, 'processed', 'raw' ) '*'] );
        if ~isempty(d)
            p = [d(1).folder filesep d(1).name];
            d = dir([p filesep '*scan.log']);
            scan_log = [d.folder filesep d.name];
        end
    end
    if ~exist( scan_log, 'file' )
        fprintf( '\nScan log not found!\n' )
        %break
    else
        %end
    
        fid = fopen( scan_log );
        cell_of_lines = textscan( fid, '%s', 'Delimiter', {'\n', '\r'} );
        cell_of_lines = cell_of_lines{1};
        fclose( fid );
        s(nn).scan_log = scan_log;
        % s_stage_z
        s_stage_z = [];
        for ll = 1:numel( cell_of_lines )
            t = regexp( cell_of_lines{ll}, 's_stage_z' );
            if t
                cc = textscan( cell_of_lines{ll}, '%s', 'Delimiter', {'=',',',' '}, 'CollectOutput', 1, 'MultipleDelimsAsOne', 1);
                s_stage_z = str2double( cc{1}{2} ) * 1000;
                break
            end
        end
    end
    if isempty(s_stage_z)
        d = dir([regexprep( full_path, 'processed', 'raw' ) '*'] );
        p = [d(1).folder filesep d(1).name];
        d = dir([p filesep '*nexus.h5']);
        h5_log = [d.folder filesep d.name];
        s_stage_z = 1000*h5read(h5_log,'/entry/scan/setup/s_stage_z');
    end
    s(nn).s_stage_z = s_stage_z;
    % Print parameters
    height = effective_pixel_size / bin * im_shape_raw(2);
    fprintf( '\n%s:', name )
    fprintf( '\n  raw_shape: %u %u', im_shape_raw )
    fprintf( '\n  raw_roi: %u %u', raw_roi )
    fprintf( '\n  bin: %u', bin )
    fprintf( '\n  effective pixelsizse: %.2f mu', effective_pixel_size )
    fprintf( '\n  s_stage_z: %.1f mu', s_stage_z )
    fprintf( '\n  full camera height (no raw_roi): %.1f mu', height)
    % Read volume
    vol_path = [full_path filesep scan_subfolder filesep reco_subfolder];
    im_struct = dir( [ vol_path filesep '*.tif' ] );
    im = imread( [vol_path filesep im_struct(1).name] );
    vol = zeros( [size(im) numel( im_struct )], 'single' );
    fprintf( '\n  Reading volume:' )
    num_slices = numel( im_struct );
    % Invert or
    filename_cell = {im_struct(end:-1:1).name};
    t = toc;
    if ~testing
        parfor kk = 1:num_slices
            filename = [vol_path filesep filename_cell{kk} ];
            vol(:,:,kk) = imread( filename, 'tif' );
        end
    end
    fprintf( ' done in %.0f s = %.1f min', toc - t, (toc - t) / 60 )
    s(nn).vol = vol;
    s(nn).size = size( vol );
    fprintf( '\n size: %u %u %u', s(nn).size )
end
%% Scan order: z-axis pointing top2bottom !! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if s(1).s_stage_z > s(2).s_stage_z
    top2bottom = 1;
    dirstr = 'top 2 bottom';
    stitch_order = num_scans:-1:1;
else
    top2bottom = 0;
    dirstr = 'bottom 2 top';
    stitch_order = 1:num_scans;
end
fprintf( '\nSample movement direction: %s', dirstr );
s(1).top2bottom = top2bottom;
s(1).stitch_order = stitch_order;


%% Normalize overlap
fprintf( '\nNormalizing volumes');
if normalize_overlap
    z = normalize_overlap;
    if top2bottom
        for nn = 1:num_scans -1
           m1 = mean3(s(nn).vol(:,:,1:z));
           m2 = mean3(s(nn+1).vol(:,:,end-z:end));
           for mm = 1:nn
               s(mm).vol = (m1 + m2) / (2 * m1) * s(mm).vol;
           end
           s(nn+1).vol = (m1 + m2) / (2 * m2) * s(nn+1).vol;
        end
    end
end

%% Absolute positions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf( '\nAbsolute positions: ' )
for nn = 1:num_scans
    
    bin = s(nn).bin;
    effective_pixel_size_binned = s(nn).effective_pixel_size_binned;
    
    % First absolute vertical position of image edge
    zval_first = s(nn).s_stage_z + double( s(nn).raw_roi(1) - 1) * s(nn).effective_pixel_size;
    s(nn).zval_first = zval_first;
    
    % Second absolute vertical position of lower image edge
    zval_last = s(nn).s_stage_z + double( s(nn).raw_roi(2) - 1 - bin ) * s(nn).effective_pixel_size;
    s(nn).zval_last = zval_last;
    
    % Absolute positions
    v1z = size( s(nn).vol, 3 );
    zval = zval_first + ( 0:v1z - 1 ) * effective_pixel_size_binned;
    s(nn).zval = zval;
    
    % Print info
    fprintf( '\n %2u. volume. first val: %.3f (%.3f) mm = %i (%i) voxel', nn, zval(1)/1000, zval_first/1000,  round(zval(1)/effective_pixel_size_binned), round(zval_first /effective_pixel_size_binned ));
    fprintf( '\n %2u. volume.  last val: %.3f (%.3f) mm = %i (%i) voxel', nn, zval(end)/1000, zval_last/1000, round(zval(end)/effective_pixel_size_binned), round(zval_last/effective_pixel_size_binned )  );
    fprintf( ', diff: %.3f (%.3f) mm', (zval(end) - zval(1))/1000, (zval_last - zval_first)/1000)
    if nn > 1 && nn <= num_scans
        if top2bottom
            if s(nn-1).zval_first > s(nn).zval_last
                cprintf('Red', '\n No overlap')
                keyboard
            else
                overlap = s(nn).zval_last - s(nn-1).zval_first;
                overlap_pixel = floor(overlap /effective_pixel_size_binned);
                s(nn-1).overlap2.mu = overlap;
                s(nn-1).overlap2.pixel = overlap_pixel;
                s(nn).overlap1.mu = overlap;
                s(nn).overlap1.pixel_binned = overlap_pixel;
                fprintf('\n overlap: %f mu, %i pixel binned', overlap, overlap_pixel)
            end
        else
            if s(nn - 1).zval_last < s(nn).zval_first
                cprintf('Red', '\n No overlap')
                keyboard
            else
                overlap = s(nn - 1).zval_last - s(nn).zval_first;
                overlap_pixel = floor(overlap /effective_pixel_size_binned);
                s(nn-1).overlap2.mu = overlap;
                s(nn-1).overlap2.pixel = overlap_pixel;
                s(nn).overlap1.mu = overlap;
                s(nn).overlap1.pixel_binned = overlap_pixel;
                fprintf('\n overlap: %f mu, %i pixel binned', overlap, overlap_pixel)
            end
        end
    end
end

%% Remove offset %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
zoffset = min( [s(:).zval] );
fap = figure( 'Name', 'Absolute Positions' );
lgnd = cell( [1 num_scans] );
fprintf( '\nAbsolute positions without offset: ' )
for nn = 1:num_scans
    s(nn).zval = s(nn).zval - zoffset;
    s(nn).zval_first = s(nn).zval_first - zoffset;
    s(nn).zval_last = s(nn).zval_last - zoffset;
    
    % Print info
    fprintf( '\n %2u. volume. first val: %f (%f) | %i %i', nn, s(nn).zval(1), s(nn).zval_first,  round(s(nn).zval(1)/effective_pixel_size_binned), round(s(nn).zval_first /effective_pixel_size_binned ));
    fprintf( '\n %2u. volume.  last val: %f (%f) | %i %i', nn, s(nn).zval(end), s(nn).zval_last, round(s(nn).zval(end)/effective_pixel_size_binned), round(s(nn).zval_last/effective_pixel_size_binned )  );
    fprintf( ', diff: %f (%f)', s(nn).zval(end) - s(nn).zval(1), s(nn).zval_last - s(nn).zval_first)
    
    % Plot absolute positions
    figure( fap )
    hold on
    plot( s(nn).zval )
    lgnd{nn} = sprintf( 'volume %u. length: %u', nn, length( zval ) ) ;
end
title( sprintf( '%s scanning', dirstr ) )
legend( lgnd )
ylabel('absolute shift')
xlabel('projection number')
if ~top2bottom
    set(gca,'ydir','reverse')
end
filename = sprintf('%s_absolute_positions.png',stitched_volume_figure);
saveas(fap,filename);

%% Noise cut level & Plot vertical cuts
yy = round( size( s(1).vol, 2 ) / 2 );
xx = round( size( s(1).vol, 1 ) / 2 );
hxz = figure( 'Name', sprintf( 'xz %s', dirstr ) );
hyz = figure( 'Name', sprintf( 'yz %s', dirstr ) );
if strcmp(noisecut,'old')
    hroi = figure( 'Name', 'vertical ROI' );
end
for nn = 1:num_scans
    if nn < num_scans && strcmp(noisecut,'new')
        hnc = figure( 'Name', sprintf('noise cut for overlap region %u', nn) );
    end
    
    % XZ
    xz = ( squeeze( s(nn).vol(:,yy,:) ) );
    figure(hxz)
    subplot(1,num_scans,stitch_order(nn))
    imsc( xz )
    zz = size( s(nn).vol, 3 );
    ztcks = 1:round( zz / 10 ):zz;
    ztcks_label = round( s(nn).zval(ztcks) );
    xtickangle(90)
    xticks( ztcks )
    xticklabels( ztcks_label )
    axis equal tight
    title( sprintf( 'volume %u', nn ) )
    % YZ
    yz = ( squeeze( s(nn).vol(xx,:,:) ) );
    figure(hyz)
    subplot(1,num_scans,stitch_order(nn))
    imsc( yz )
    xticks( ztcks )
    xtickangle(90)
    xticklabels( ztcks_label )
    axis equal tight
    title( sprintf( 'volume %u', nn ) )
    switch noisecut
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        case 'none'
            fprintf('\n No noise related volume cropping is applied')
            % Set noisecut to full volume
            if top2bottom
                %% top 2 bottom %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % 1 bottom
                % 1 top
                %   2 bottom
                %   2 top
                if nn == 1
                    s(nn).noisecut2 = s(nn).size(3);
                end
                if nn >= 1 && nn < num_scans
                    ind1 = 1;
                    ind2 = s(nn+1).size(3);
                    s(nn).noisecut1 = ind1;
                    s(nn+1).noisecut2 = ind2;
                end
                if nn == num_scans
                    s(nn).noisecut1 = 1;
                end
            else
                %% bottom 2 bop %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % 1 top
                % 1 bottom
                % 2 top
                % 2 bottom
                if nn == 1
                    s(nn).noisecut1 = 1;
                end
                if nn >= 1 && nn < num_scans
                    ind1 = s(nn).size(3);
                    ind2 = 1;
                    s(nn).noisecut2 = ind1;
                    s(nn+1).noisecut1 = ind2;
                end
                if nn == num_scans
                    s(nn).noisecut2 = s(nn).size(3);
                end
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        case 'old'
            % ROI regarding noise level
            im = xz(end-100:end,:);
            l = ( std( stdfilt( im ), 1, 1) ) ;
            l = l(:);
            l = l(end:-1:1);
            % Filter
            sigma = 15;
            sz = 3*length( l );    % length of gaussFilter vector
            x = linspace(-sz / 2, sz / 2, sz);
            gaussFilter = exp(-x .^ 2 / (2 * sigma ^ 2));
            gaussFilter = gaussFilter / sum (gaussFilter); % normalize
            lf = conv( padarray( l, [length(l) 0], 'symmetric', 'both' ), gaussFilter, 'same');
            lf = lf(length(l) + (1:length(l) ) );
            
            lfx = round( length(l)/4):round( length(l)*3/4);
            t = stitch_level_fac * mean( lf(lfx) );
            % Lower index
            v1 = 1;
            while 1
                if lf(v1) > t
                    v1 = v1 + 1;
                else
                    break
                end
            end
            s(nn).noisecut1 = v1;
            %Higher index
            v2 = length( l );
            while 1
                if lf(v2) > t
                    v2 = v2 - 1;
                else
                    break
                end
            end
            s(nn).noisecut2 = v2;
            
            % Plot
            figure( hroi )
            subplot(1,num_scans,nn)
            plot( [l, lf, t*ones(size(l))] )
            title( sprintf( 'xz cut level. left: %u, right: %u', v1, v2 ) )
            axis tight
        case 'new'
            %if nn > 1 && nn <= num_scans
            
            if top2bottom
                %% top 2 bottom %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % 1 bottom
                % 1 top
                %   2 bottom
                %   2 top
                if nn == 1
                    s(nn).noisecut2 = s(nn).size(3);
                end
                if nn >= 1 && nn < num_scans
                    v1m = squeeze( mean( s(nn).vol, [1 2]));
                    v1s = squeeze(std( s(nn).vol, 0, [1 2]));
                    v2m = squeeze(mean( s(nn+1).vol, [1 2]));
                    v2s = squeeze(std( s(nn+1).vol, 0, [1 2]));
                    v1r = v1s./v1m;
                    v2r = v2s./v2m;
                    figure(hnc)
                    subplot(num_scans-1,1,nn)
                    plot( [v2r; v1r],'.')
                    ylabel( 'std/mean')
                    xlabel( 'slice' )
                    axis tight
                    overlap_pixel = floor((s(nn+1).zval_last - s(nn).zval_first ) /effective_pixel_size_binned);
                    % set crop level to areas with similiar SNR
                    ind1 = 1;
                    ind2 = s(nn+1).size(3);
                    %increase/decrease indices if SNR is larger than the
                    %other
                    for ll = 1:overlap_pixel - 1
                        if v1r(ind1) >= v2r(ind2)
                            ind1 = ind1 + 1;
                        else
                            ind2 = ind2 - 1;
                        end
                    end
                    % Set noise cuts
                    s(nn).noisecut1 = ind1;
                    s(nn+1).noisecut2 = ind2;
                    % Plot crop points
                    hold on
                    plot( ind2, v2r(ind2), 'ro')
                    plot( s(2).size(3) + ind1, v1r(ind1), 'ro')
                    filename = sprintf( '%s_figure_noisecut_overlap%02u.png', stitched_volume_figure, nn);
                    saveas( hnc,  filename);
                end
                if nn == num_scans
                    s(nn).noisecut1 = 1;
                end
            else
                %% bottom 2 bop %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % 1 top
                % 1 bottom
                % 2 top
                % 2 bottom
                if nn == 1
                    s(nn).noisecut1 = 1;
                end
                if nn >= 1 && nn < num_scans
                    v1m = squeeze( mean( s(nn).vol, [1 2]));
                    v1s = squeeze(std( s(nn).vol, 0, [1 2]));
                    v2m = squeeze(mean( s(nn+1).vol, [1 2]));
                    v2s = squeeze(std( s(nn+1).vol, 0, [1 2]));
                    v1r = v1s./v1m;
                    v2r = v2s./v2m;
                    figure(hnc)
                    subplot(num_scans-1,1,nn)
                    plot( [v1r; v2r],'.')
                    ylabel( 'std/mean')
                    xlabel( 'slice' )
                    axis tight
                    overlap_pixel = floor(( s(nn).zval_last - s(nn+1).zval_first) /effective_pixel_size_binned);
                    % set crop level to areas with similiar SNR
                    ind1 = s(nn).size(3);
                    ind2 = 1;
                    %increase/decrease indices if SNR is larger than the
                    %other
                    for ll = 1:overlap_pixel - 1
                        if v1r(ind1) >= v2r(ind2)
                            ind1 = ind1 - 1;
                        else
                            ind2 = ind2 + 1;
                        end
                    end
                    % Set noise cuts
                    s(nn).noisecut2 = ind1;
                    s(nn+1).noisecut1 = ind2;
                    % Plot crop points
                    hold on
                    plot( ind1, v1r(ind1), 'ro')
                    plot( s(1).size(3) + ind2, v2r(ind2), 'ro')
                    filename = sprintf( '%s_figure_noisecut_overlap%02u.png', stitched_volume_figure, nn);
                    saveas( hnc,  filename);
                end
                if nn == num_scans
                    s(nn).noisecut2 = s(nn).size(3);
                end
            end
    end
end
drawnow
fprintf('\n')

%% Overlap %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for nn = 1:num_scans - 1
    % Absolute z ranges of subsequent volumes
    zval1 = s(nn).zval;
    zval2 = s(nn+1).zval;
    % Find overlap (within noisecut ROI if possible)
    if top2bottom
        % First / upper overlap index and value in first volume
        z1ol_pos = s(nn).noisecut1;
        z1ol_val = zval1(z1ol_pos);
        % Check if first/upper overlap edge of first volume is below lower edge of second volume
        if z1ol_val > zval2(end)
            z1ol_pos = 1;%zval1(1);
            z1ol_val = zval1(z1ol_pos);
        end
        % Second / lower overlap index and value in second volume
        z2ol_pos = s(nn+1).noisecut2;
        z2ol_val = zval2(z2ol_pos);
        % Check if second/lower overlap edge of second volume is above first
        % volume
        if z2ol_val < zval1(1)
            %z2ol_pos = zval2(end);
            %z2ol_val = zval2(z2ol_pos);
            z2ol_pos = numel(zval2);
            z2ol_val = zval2(z2ol_pos);
        end
    else
        % Lower overlap edge in first volume
        z1ol_pos = s(nn).noisecut2;
        z1ol_val = zval1(z1ol_pos);
        % Check if lower overlap edge of first volume is above upper edge
        % of second volume
        if z1ol_val < zval2(1)
            z1ol_pos = numel(zval1);%(end);
            z1ol_val = zval1(z1ol_pos);
        end
        % Upper overlap edge in second volume
        z2ol_pos = s(nn+1).noisecut1;
        z2ol_val = zval2(z2ol_pos);
        % Check if the lower end of the ROI of the second volume is still
        % larger than the upper end of the first volume
        if z2ol_val > zval1(end)
            %z2ol_pos = zval2(1);
            %z2ol_val = zval2(z2ol_pos);
            z2ol_pos = 1;
            z2ol_val = zval2(z2ol_pos);
        end
    end
    % Stitch level
    zstitch = ( z1ol_val + z2ol_val ) / 2;
    [~, z1stitch_pos] = min( abs( zval1 - zstitch ) );
    [~, z2stitch_pos] = min( abs( zval2 - zstitch ) );
    z1stitch_val = zval1( z1stitch_pos );
    z2stitch_val = zval2( z2stitch_pos );
    if top2bottom
        
        % Check monotonicity
        while z1stitch_val < z2stitch_val
            z1stitch_pos = z1stitch_pos + 1;
            z1stitch_val = zval1( z1stitch_pos );
        end
        
        % Stich indices
        s(nn).zstitch_pos1 = z1stitch_pos;
        if nn == num_scans - 1
            s(nn).zstitch_pos2 = size( s(nn).vol, 3 );
        end
        if nn == 1
            s(nn+1).zstitch_pos1 = 1;
        end
        s(nn+1).zstitch_pos2 = z2stitch_pos;
    else
        % Check monotonicity
        while z1stitch_val > z2stitch_val
            z2stitch_pos = z2stitch_pos + 1;
            z2stitch_val = zval2( z2stitch_pos );
        end
        % Stich indices
        if nn == 1
            s(nn).zstitch_pos1 = 1;
        end
        s(nn).zstitch_pos2 = z1stitch_pos;
        s(nn+1).zstitch_pos1 = z2stitch_pos;
        if nn == num_scans - 1
            s(nn+1).zstitch_pos2 = size( s(nn+1).vol, 3 );
        end
    end
    % Print overlap info
    fprintf( '\n %u. Overlap:', nn )
    if top2bottom
        fprintf( '\n   %-30s %f', '1. vol lower edge', zval1(end) )
        fprintf( '\n   %-30s %f', '2. vol lower edge',  zval2(end) )
        fprintf( '\n   %-30s %f', '2. vol lower noisecut2', zval2( s(nn+1).noisecut2 ) )
        fprintf( '\n   %-30s %f', 'stitch level', zstitch )
        fprintf( '\n   %-30s %f', '1. vol upper noisecut1', zval1( s(nn).noisecut1 ) )
        fprintf( '\n   %-30s %f', '1. vol upper edge', zval1(1) )
        fprintf( '\n   %-30s %f', '2. vol upper edge', zval2(1) )
    else
        fprintf( '\n   %-30s %f', '1. vol upper edge', zval1(1) )
        fprintf( '\n   %-30s %f', '2. vol upper edge', zval2(1) )
        fprintf( '\n   %-30s %f', '2. vol noisecut1', zval2( s(nn+1).noisecut1 ) )
        fprintf( '\n   %-30s %f', 'stitch level', zstitch )
        fprintf( '\n   %-30s %f', '1. vol lower noisecut2', zval1( s(nn).noisecut2) )
        fprintf( '\n   %-30s %f', '1. vol lower edge', zval1(end) )
        fprintf( '\n   %-30s %f', '2. vol lower edge', zval2(end) )
    end
    % Print stitch level info
    fprintf( '\n Stitch level:')
    fprintf( '\n   %-30s value: %f', 'stitch level.', zstitch )
    fprintf( '\n   %-30s value: %f, index: %u', 'stitch level volume 1.', z1stitch_val, z1stitch_pos )
    fprintf( '\n   %-30s value: %f, index: %u', 'stitch level volume 2.', z2stitch_val, z2stitch_pos )
    fprintf( '\n')
    % Stitch values, positions and ranges
    if top2bottom
        zz1 = z1stitch_pos:length(zval1);
        zz2= 1:z2stitch_pos;
    else
        zz1 = 1:z1stitch_pos;
        zz2 = z2stitch_pos:length( zval2 );
    end
    % Plot
    h = figure('Name', sprintf( 'vertically stitched volume part %u', nn ) );
    % XZ
    imx1 = squeeze( s(nn).vol(:,yy,zz1) );
    imx2 = squeeze( s(nn+1).vol(:,yy,zz2) );
    % YZ
    imy1 = squeeze( s(nn).vol(xx,:,zz1) );
    imy2 = squeeze( s(nn+1).vol(xx,:,zz2) );
    if top2bottom
        xzstitch = rot90( cat( 2, imx2, imx1 ), -1 );
        yzstitch = rot90( cat( 2, imy2, imy1 ), -1 );
    else
        xzstitch = rot90( cat( 2, imx1, imx2 ), -1 );
        yzstitch = rot90( cat( 2, imy1, imy2 ), -1 );
    end
    subplot(1,2,1)
    imsc( FilterHisto(xzstitch ))
    axis equal tight
    title( 'xz plane' )
    subplot(1,2,2)
    imsc( FilterHisto(yzstitch ))
    axis equal tight
    title( 'yz plane' )
    drawnow
    % Save figures and ortho views
    filename = sprintf( '%s_figure_part%02u.png', stitched_volume_figure, nn);
    saveas( h,  filename);
    filename = sprintf( '%s_xz_part%02u.tif', stitched_volume_figure, nn);
    write32bitTIFfromSingle( filename, xzstitch );
    filename = sprintf( '%s_yz_part%02u.tif', stitched_volume_figure, nn);
    write32bitTIFfromSingle( filename, yzstitch );
end
%% z index range for stitching
fprintf( '\nz indices for stitching:' )
for nn = 1:num_scans
    
    z1 = s(nn).zstitch_pos1;
    if isempty(z1)
        z1 = 1;
    end
    z2 = s(nn).zstitch_pos2;
    if isempty(z2)
        z2 = size(s(nn).vol,3);
    end
    s(nn).zstitch_range = z1:z2;
    fprintf( '\n %2u. vol: %5u, %5u (of %5u)', nn, z1, z2, s(nn).size(3) )
end

%% Save stitched volume
CheckAndMakePath(stitched_volume_path)
% if ~iscell( scan_path2 )
%     [~,name]=fileparts(scan_path2);
% else
%     [~,name]=fileparts(scan_path{end});
% end
if top2bottom
    index_offset = 0;
    for nn = 1:num_scans
        zrange = s(nn).zstitch_range; % slices taken
        num_slices = numel(zrange);
        zmax = zrange(end); % highest slice number used
        vol = s(nn).vol;
        t = toc;
        fprintf( '\nSaving stitched volume slices:' )
        if ~testing
            parfor ll = zrange % loop over slices taken
                % create properly increasing file index
                ind =  1 + (zmax - ll) + index_offset
                filename = sprintf( '%s/%s_%06u.tif', stitched_volume_path, stitched_name, ind);
                im = vol(:,:,ll);
                write32bitTIFfromSingle( filename, im )
            end
        end
        fprintf( ' done in %.0f s = %.1f min', toc - t, (toc - t) / 60 )
        index_offset = index_offset + num_slices; % needed for image index
    end
else
    index_offset = 0;
    for nn = 1:num_scans
        zrange = s(nn).zstitch_range;
        num_slices = numel(zrange);
        zmin = zrange(1);
        vol = s(nn).vol;
        fprintf( '\nSaving stitched volume slices:' )
        t = toc;
        if ~testing
            parfor ll = zrange
                ind =  1 + (ll - zmin) + index_offset
                filename = sprintf( '%s/%s_%06u.tif', stitched_volume_path, stitched_name, ind);
                im = vol(:,:,ll);
                write32bitTIFfromSingle( filename, im )
            end
        end
        fprintf( ' done in %.0f s = %.1f min', toc - t, (toc - t) / 60 )
        index_offset = index_offset + num_slices;
    end
end
fprintf( '\nSTITCHING FINISHED in %.0f s = %.1f min\n\n', toc, toc / 60 )
