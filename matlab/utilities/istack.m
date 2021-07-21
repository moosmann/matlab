function istack( str, transpose_plots )
% Display ortho slices from tomo sequence that are created automatically
% during reconstruction.

if nargin < 1
    folder = pwd;
    [p,f] = fileparts( folder );
    str = [p filesep f(1:end-1) '*'];
end
if nargin < 2
    transpose_plots = 0;
end
if nargin < 3
    num_im = 4;
end
if nargin < 4
    skip = 10;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist( 'f', 'var' )
        [~,f] = fileparts( str );
end
fprintf( '\nscan: %s', f )
fprintf( '\nstring: %s', str )
s = dir( str );

steps = numel( s );
fprintf( '\n steps found: %u', steps )

if steps > num_im
    ind = unique( [ round( 1:steps/num_im:steps), steps ] );
    s = s(ind);
    fprintf( '\n indices used: ' )
    fprintf( ' %u', ind )
end

figure('Name', f ,  'WindowState', 'maximized');
c = 0;
h = numel(s);
for n = 1:h
    
    imx = imread( sprintf( '%s/%s/reco/reco_xMid.tif', s(n).folder, s(n).name ) );
    if transpose_plots
        subplot( 3, h, n)
    else
       subplot( h, 3, 1 + c)
    end
    imsc(flipud( imx(1:skip:end,1:skip:end)) )
    title( ['x ' s(n).name], 'Interpreter', 'none' )
    axis tight equal
    
    imy = imread( sprintf( '%s/%s/reco/reco_yMid.tif', s(n).folder, s(n).name ) );
    if transpose_plots
        subplot( 3, h, n + h)
    else
       subplot( h, 3, 2 + c)
    end
    imsc( flipud( imy(1:skip:end,1:skip:end)) )
    title( ['y ' s(n).name], 'Interpreter', 'none' )
    axis tight equal
    
    imz = imread( sprintf( '%s/%s/reco/reco_zMid.tif', s(n).folder, s(n).name ) );
    if transpose_plots
        subplot( 3, h, n + 2 * h)
    else
       subplot( h, 3, 3 + c)
    end
    imsc( imz(1:skip:end,1:skip:end))
    title( ['z ' s(n).name], 'Interpreter', 'none' )
    axis tight equal
    
    if transpose_plots
         
    else
        c = c + 3;
    end
end

fprintf( '\n' )

