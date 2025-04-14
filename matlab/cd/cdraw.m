function cdraw

fid = fopen( sprintf( '%s/path_to_raw', userpath ), 'r');
cd(fscanf( fid,'%s'))
fclose( fid );