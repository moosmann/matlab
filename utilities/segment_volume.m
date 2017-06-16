function [vol, out] = segment_volume( vol, nbins, visual_output, verbose )

if nargin < 2
    nbins = 2^10;
end
if nargin < 3
   visual_output = 1;
end
if nargin < 4
    verbose = 1;
end

%% Main %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t = toc;
PrintVerbose(verbose, 'Segmentation:')    
    
d3 = size( vol, 3);
d3h = round( d3 / 2 );

% Preallocation
hs = zeros(nbins,3);
m1 = zeros(d3,1);
m2 = m1;
m3 = m2;
mm = m1;
st = m1;
max1_val = m1;
max2_val = m1;
max1_pos = m1;
max2_pos = m1;
fwhm = m1;
t1_val = m1;
t1_pos = m1;
x1 = m1;

x = 1:nbins;

% mask for central disc (ellipse)
[~, m] = MaskingDisc( vol(:,:,d3h), 0.95, 0);

vol_min = min3( bsxfun(@times, vol, m) );
vol_max = max3( bsxfun(@times, vol, m) );

% Edges of bin interval
edges = (0:nbins) / nbins * (vol_max - vol_min) + vol_min;

parfor nn = 1:d3
    s = vol(:,:,nn);
    h = histcounts( s, edges);
    hs(:,nn) = h;
    
    % Maximum: position and value
    [max1_val(nn), max1_pos(nn)] = max( h );
    
    % FWHM
    [~, p] = min( abs( h - max1_val(nn)/2 ) );
    fwhm(nn) = abs( max1_pos(nn) - p );
    sig = fwhm(nn) / ( 2 * sqrt( 2 * log(2) ) );
    
    % Search 2nd maximum
    x1(nn) = ceil( max1_pos(nn) + 11*sig);
    [max2_val(nn), max2_pos(nn)] = max( h(x1(nn):end) );
    max2_pos(nn) = max2_pos(nn) + x1(nn) - 1 ;
        
    % Momenta
    m1(nn) = sum( x .* h ) / sum(h);
    m2(nn) = (sum( x.^2 .* h ) / sum(h) )^(1/2);
    m3(nn) = (sum( x.^3 .* h ) / sum(h) )^(1/3);
    
    % Mean and standard deviations
    mm(nn) = mean2( s );
    st(nn) = std2( s );
    
    % Threshold
    t1_pos(nn) = round( (max1_pos(nn) + max2_pos(nn))/2 );
    t1 = sum( edges( t1_pos(nn) + (0:1) ) )/2;
    t1_val(nn) = t1;
    
    % segmentation: background = 0, screw = 1
    vs = zeros( size( s ) , 'uint8');
    vs(s >= t1) = 1;
    
    %% Try simple segmentaiton of Gd dots
    [h2,edges2] = histcounts( s(logical(vs)) );
    [~,max_pos] = max( h2 );
    t2 = edges2(floor(1.8*max_pos));
    vs( s >= t2 ) = 2;
        
    %s(s > t1) = 1;        
    vol(:,:,nn) = vs;
    
end

%% output struct
out.histo = hs;
out.edges = edges;
out.m1 = m1;
out.m2 = m2;
out.mean = mm;
out.std = st;
out.max.val = max1_val;
out.max.pos = max1_pos;
out.fwhm = fwhm;
out.max2.val = max2_val;
out.max2.pos = max2_pos;
out.t1.val = t1_val;
out.t1.pos = t1_pos;
out.x1 = x1;

%% Show plots
if visual_output
    %figure('Name', 'volume histogram');
    m = 2;
    n = 4;
    subplot(m,n,1)
    imsc( log(1 + ( hs ) ) );
    axis equal tight
    title(sprintf('histo stack'))    
    set(gca,'YDir','normal')
    
    subplot(m,n,2)
    Y = [m1, m2, m3, max1_pos];
    plot( Y )
    axis tight
    legend( '1st momentum', '2nd momentum', '3rd momentum', 'maximum position' )   
    title(sprintf('histogram: momentum & maximum'))    
     
    subplot(m,n,3)
    Y = [mm, st];
    plot( Y )
    axis tight
    legend( 'mean', 'std')    
    title(sprintf('slices: mean & std'))    
    
    subplot(m,n,4)
    Y = [max1_val, max2_val];
    plot( Y )
    axis tight
    legend( 'maximum 1', 'maximum 2')
    title(sprintf('histogram: maxima values'))    
    
    subplot(m,n,5)
    Y = [max1_pos, max2_pos, x1];
    plot( Y )
    axis tight
    legend( 'maximum 1', 'maximum 2', 'hist cut')
    title( 'histogram: maxima position')   
    
    subplot(m,n,6)
    Y = [t1_pos];
    plot( Y )
    axis tight
    legend( 'threshold')    
    title( 'threshold position')
        
    subplot(m,n,7)
    Y = [t1_val];
    plot( Y )
    axis tight
    legend( 'threshold value')    
    title( 'threshold value')
    
    subplot(m,n,8)
    imsc( vol(:,:,d3h))    
    txt = text(20,40, sprintf( '%u', d3h) );
    txt.Color = [1 1 1];
    txt.FontSize = 16;
    axis image    
    title( 'segmentaion')
    
    drawnow
end

PrintVerbose(verbose, ' done in %.0f s (%.2f min).', toc - t , (toc - t) / 60)



