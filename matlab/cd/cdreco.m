function cdreco

fid = fopen( sprintf( '%s/path_reco', userpath ), 'r');
cd( fscanf( fid, '%s' ) )
fclose( fid );
