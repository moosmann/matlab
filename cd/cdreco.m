function cdreco

fid = fopen( sprintf( '%s/experiments/p05/path_to_latest_reco', userpath ), 'r');
cd( fscanf( fid, '%s' ) )
fclose( fid );
