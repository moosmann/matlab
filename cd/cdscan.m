function cdscan

fid = fopen( sprintf( '%s/experiments/p05/path_to_latest_scan', userpath ), 'r');
cd( fscanf( fid, '%s' ) )
fclose( fid );