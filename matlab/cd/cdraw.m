function cdraw

fid = fopen( sprintf( '~/path_to_raw'),'r');
cd(fscanf( fid,'%s'))
fclose( fid );