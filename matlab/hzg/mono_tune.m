fn =  '/asap3/petra3/gpfs/common/p05/jm/matlab/mono_tune.txt';
fid = fopen( fn );
c = textscan( fid, '%*s%f\n%*s%f\n%*s=%f%*s=%f%*s=%f');%, 'Delimiter', {'\n', '\r'} )
fclose( fid );

b1 = c{1};
b2 = c{2};
p2 = c{3};
mi = c{4};
fwhm = c{5};
%plot3(b1,b2,mi,'o')
%scatter3(b1,b2,mi)
z = 100 * normat( mi) + 1;
c = normat(fwhm);
figure( 'Name', 'scatter' )
 scatter(b1,b2,z,c,'filled')
 
figure( 'Name', 'scatter 2nd bend 2nd pitch' )
 scatter(b2,p2,z,c,'filled')
xlabel( '2nd bend' )
ylabel( '2nd pitch' )
 
figure( 'Name', 'scatter3' )
scatter3(b1,b2,p2,z,c,'filled')
colormap copper
xlabel( '1st bend' )
ylabel( '2nd bend' )
zlabel( ' 2nd pitch' )
%pos_p07_dcm_1st_bend     0.019997636
%pos_p07_dcm_2nd_bend     -0.48598097
%2nd_pitch =        2.1439649 MaxInt =        13053.769  FWHM =     0.0055839391
