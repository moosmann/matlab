% Create multi tiff from load sequence
% Optionally crop and compress

scan_name = 'syn13_55L_Mg10Gd_12w_load';
par_path = '/asap3/petra3/gpfs/p05/2017/data/11003950/processed';
CheckTrailingSlash( par_path )

outpath = [ par_path(1:end-4) '/processed/' scan_name '/'];
CheckTrailingSlash( outpath );

fprintf( '\nscan_path : %s', par_path)
fprintf( '\noutpath : %s', outpath)
 
sequ_indices = [0];

for nn = sequ_indices
    
    scan_path = sprintf( '%s%s_%02u', par_path, scan_name, nn);        
    fprintf( '\n%s', scan_path )
            
end

fprintf( '\n' )