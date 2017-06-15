function cdraw

fid = fopen( sprintf( '%s/experiments/p05/path_to_latest_raw', userpath ), 'r');
cd( fscanf( fid, '%s' ) )
fclose( fid );