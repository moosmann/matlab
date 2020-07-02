function cdsavedworkingdirectory

fid = fopen( '~/.savedworkingdirectory', 'r');
swd = fscanf( fid, '%s' );
cd( swd )
fclose( fid );
