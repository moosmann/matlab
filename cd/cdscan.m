function cdscan

fid = fopen( sprintf( '%s/experiments/p05/pathtolastscan', userpath ), 'r');
cd( fscanf( fid, '%s' ) )
fclose( fid );