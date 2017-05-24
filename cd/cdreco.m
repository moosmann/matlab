function cdreco

fid = fopen( sprintf( '%s/experiments/p05/pathtolastreco', userpath ), 'r');
cd( fscanf( fid, '%s' ) )
fclose( fid );
