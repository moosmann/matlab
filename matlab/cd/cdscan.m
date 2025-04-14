function cdscan

fid = fopen( sprintf( '%s/path_to_scan', userpath ), 'r');
cd( fscanf( fid, '%s' ) )
fclose( fid );