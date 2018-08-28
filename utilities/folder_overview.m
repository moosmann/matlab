function folder_overview( pattern )

if nargin < 1
    pattern = './*';
end

%% Main

d = dir( pattern );
isub = [d(:).isdir];
d = {d(isub).name};
while d{1}(1) == '.'
    d(1) = [];
end

for nn = 1:numel( d )
   
    folder = d{nn};
    [~,num_files] = unix( sprintf( 'ls -lA %s | wc -l;', folder ) );
    
    num_files = str2double( num_files );
    %num_files = numel( dir( folder ) ) - 2;
    
    fprintf( '\n%3u %6u %s', nn, num_files, folder )
end

fprintf('\n')


