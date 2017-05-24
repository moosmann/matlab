function cdraw

fid = fopen( sprintf( '%s/experiments/p05/pathtolastraw', userpath ), 'r');
cd( fscanf( fid, '%s' ) )
fclose( fid );