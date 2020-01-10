% Magnification vs distance using KIT camera

% vertical MTF, horizontal edge
mtf_vert = [ 
80   4.88725
200  4.89246
400  4.9040126
600  4.9130959
800  4.9237997
1000 4.9322013
1200 4.9425318
];

% horizontal MTF, vertical edge
mtf_hor = [
1200 4.9375306
1000 4.9281522
 800 4.9184230
 600 4.9091415
 400 4.8993784
 200 4.8899764
  80 4.8843443
];

% Plot magn
x = mtf_vert(:,1);
mv = mtf_vert(:,2);
mh = mtf_hor(end:-1:1,2);
figure( 'Name', 'magnification vs distance' )
plot( x, mv, x, mh)
xlabel( ' z / mm' )
ylabel( ' magn ' )
legend( {'vertical' 'horizontal'} )

% Fit 
fitType = 'poly1';
fov = fit( x, mv, fitType);
foh = fit( x, mh, fitType);

% Print fit parameters
disp( fov )
disp( foh )

% extrapolate source distance
fprintf( '\nsource distance extrapolation : \n hor : %g m\n vert : %g m', foh.p2 / foh.p1 / 1000, fov.p2 / fov.p1 / 1000)


 mh = @(x) foh.p1 * x + foh.p2;
 mv = @(x) fov.p1 * x + fov.p2;
 
 % Data reco
 scan_path = '/asap3/petra3/gpfs/p07/2019/data/11007454/processed/bmc06_tooth1';
 % empty slices
 [3996 3997 4062 4063 ]
 %noisy slices, probably beam dump / missing information / etc
 4516:5259 % approx
 
