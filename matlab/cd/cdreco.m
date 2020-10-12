function cdreco

fid = fopen( sprintf( '%s/path_to_reco', userpath ), 'r');
cd( fscanf( fid, '%s' ) )
fclose( fid );
