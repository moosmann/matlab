function istack( str )

if nargin < 1
    folder = pwd;
    [p,f] = fileparts( folder );
    str = [p filesep f(1:end-1) '*'];
end

if ~exist( 'f', 'var' )
        [~,f] = fileparts( str );
end
fprintf( '\nscan: %s', f )
fprintf( '\nstring: %s', str )
s = dir( str );

figure('Name', f ,  'WindowState', 'maximized');
c = 0;
h = numel(s);
for n = 1:h
    
    imx = imread( sprintf( '%s/%s/reco/reco_xMid.tif', s(n).folder, s(n).name ) );
    subplot( h, 3, 1 + c)
    imsc(flipud( imx) )
    title( ['x ' s(n).name], 'Interpreter', 'none' )
    axis tight equal
    
    
    imy = imread( sprintf( '%s/%s/reco/reco_yMid.tif', s(n).folder, s(n).name ) );
    subplot( h, 3, 2 + c)
    imsc( flipud( imy) )
    title( ['y ' s(n).name], 'Interpreter', 'none' )
    axis tight equal
    
    imz = imread( sprintf( '%s/%s/reco/reco_zMid.tif', s(n).folder, s(n).name ) );
    subplot( h, 3, 3 + c)
    imsc(imz)
    title( ['z ' s(n).name], 'Interpreter', 'none' )
    axis tight equal
    
    c = c + 3;
end

fprintf( '\n' )

