fprintf( '\nReading data')
tic
%p = '/asap3/petra3/gpfs/p07/2021/data/11013549/raw/hereon_114_edge_500/';
p = '/asap3/petra3/gpfs/p07/2021/data/11013549/raw/hereon_114hb_edge_1050/';

curstedgex = read_dat([p 'curstedgex.dat']);
curstedgez = read_dat([p 'curstedgez.dat']);
curstref = read_dat([p 'curstref.dat']);

imgstdark = read_dat([p 'imgstdark.dat']);
imgstedgex = read_dat([p 'imgstedgex.dat']);
imgstedgez = read_dat([p 'imgstedgez.dat']);
imgstref = read_dat([p 'imgstref.dat']);

xstagestedge= read_dat([p 'xstagestedge.dat']);
zstagestedge = read_dat([p 'zstagestedge.dat']);
toc

figure('Name', 'beam current')
plot([curstedgex; curstedgez; curstref]')
legend( {'edge x', 'edge z', 'ref'} )
axis tight

figure('Name', 'egde position')
plot([xstagestedge; zstagestedge]')
legend( {'edge x', 'edge z'} )
axis tight

%%
fprintf('\n flat field correction')
tic
bin = 5;
par.threshold_hot = 0.05;
par.threshold_dark = 0.05;
par.medfilt_neighboorhood = [5 5];
par.filter_dead_pixel = 1;
par.verbose = 0;
par.use_gpu = 1;
for n = numel( curstref ):-1:1
    dark = Binning( FilterPixel( imgstdark(:,:,n), par), bin );
    ref = Binning( FilterPixel( imgstref(:,:,n), par), bin );
    imx(:,:,n) = Binning( FilterPixel( imgstedgex(:,:,n), par), bin );
    imz(:,:,n) = Binning( FilterPixel( imgstedgez(:,:,n), par), bin );
    edgex(:,:,n) = ( imx(:,:,n) - dark ) ./ (ref - dark );
    edgez(:,:,n) = ( imz(:,:,n) - dark) ./ (ref - dark);
    
end
toc
%%
indx = 10:19;
indz = 33:49;
edgex = edgex(:,:,indx);
edgez = edgez(:,:,indz);






