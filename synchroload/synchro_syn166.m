clear p s f t l t0
scan = 'syn166_104R_Mg10Gd_4w';
p = '/asap3/petra3/gpfs/p05/2017/data/11003288/raw/';

outpath = [ p(1:end-4) '/processed/' scan '/'];
CheckTrailingSlash( p );

fprintf( 'outpath : %s', outpath)


for nn = 9:-1:0
    
    s = sprintf( '%s_%03u', scan, nn);        
    filename = sprintf( '%s%s/%stime_adc1_p05pusher_tomo.dat', p, s, s );
    fprintf( '\n%u  %s', exist(filename, 'file') ,filename )
    
    f = read_dat( filename );
    d(nn+1,:,:) = f(:,2:end-1);
        
end

calib = 4.81879;
t0 = squeeze( d(1,1,1) );
for nn = 10:-1:1
    t(:,nn) = squeeze( d(nn,1,:) );
    l(:,nn) = calib * squeeze( d(nn,2,:) ) ;    
end

%% FIGURE
fontSize = 20;
prnt =  @(filename) print(filename,'-dpng','-r300');
fgr = @(name) figure( 'Name', name, 'PaperPositionMode', 'Manual', 'PaperPosition', [0 0 32 16] );
if exist('fig1','var') && isvalid( fig1 )
    fig1 = figure( fig1 );
else
    fig1 = fgr( 'step 2');
end
clf( fig1 );
ax = axes( 'Parent', fig1, 'FontSize',fontSize );
box( ax, 'on');
hold( ax, 'all');
h = plot( l(:,2), 'Parent', ax);
axis tight
set(h, 'LineWidth',2);
xlabel('number of images','FontSize',fontSize)
ylabel('force / N','FontSize',fontSize)
filename = sprintf( '%s%s_load_step2.png', outpath, scan(1:6) );
set(gcf,'PaperPositionMode','auto')
prnt( filename )

%% FIGURE

if exist('fig2', 'var') && isvalid( fig2 )
    fig2 = figure( fig2 );
else
    fig2 = fgr( 'all steps' );
end
clf( fig2 );
ax = axes( 'Parent', fig2, 'FontSize',fontSize );
box( ax, 'on');
hold( ax, 'all');
h = plot( l(:), '.', 'Parent', ax);
axis tight
set(h, 'LineWidth',4);
xlabel('number of images','FontSize',fontSize)
ylabel('force / N','FontSize',fontSize)
filename = sprintf( '%s%s_load_stepsAll.png', outpath, scan(1:6) );
set(gcf,'PaperPositionMode','auto')
prnt( filename )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf( '\n' );